import os

REF=f"{config['species']}.{config['genome']}"

rule directories:
    output: directory(f"data/ensembl/{REF}.dna.primary_assembly.karyotypic/")
    shell: "mkdir -p data/ensembl/{REF}.dna.primary_assembly.karyotypic/"

rule hisat_genome:
    '''Build genome index for hisat2'''
    input:
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        gtf=f"data/ensembl/{REF}.{config['release']}.gff3",
    threads: 12
    output:
        idx=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.1.ht2",
        finished=f"data/ensembl/done_building_hisat_genome{REF}.txt",
    benchmark: f"data/ensembl/{REF}.hisatbuild.benchmark"
    log: f"data/ensembl/{REF}.hisatbuild.log"
    shell:
        "(hisat2-build -p {threads} data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa"
        " data/ensembl/{REF}.dna.primary_assembly.karyotypic && touch {output.finished}) &> {log}"

rule hisat2_splice_sites:
    '''Fetch the splice sites from the gene model for hisat2'''
    input: f"data/ensembl/{REF}.{config['release']}.gff3"
    output: f"data/ensembl/{REF}.{config['release']}.splicesites.txt"
    shell: "hisat2_extract_splice_sites.py {input} > {output}"

def check_sra():
    '''Check if SRAs should be downloaded'''
    docheck = 'sra' in config and config["sra"] is not None and len(config["sra"]) > 0
    return docheck

def sra_output(wildcards):
    prefix = f"{wildcards.dir}/{wildcards.sra}"
    spotfile = f"{prefix}.spot"
    fq1 = f"{prefix}_1.fastq"
    fq2 = f"{prefix}_2.fastq"
    isSingleEnd = sum([1 for line in open(spotfile)]) == 4
    return [fq1] if isSingleEnd else [fq1, fq2]

def fastqs(wildcards, returnCompressed):
    prefix = f"{wildcards.dir}/{wildcards.fq}"
    fq1gz = f"{prefix}_1.fastq.gz"
    fq2gz = f"{prefix}_2.fastq.gz"
    fq1 = fq1gz.rstrip('.gz')
    fq2 = fq2gz.rstrip('.gz')
    isSingleEnd = not os.path.exists(fq2gz)
    if returnCompressed:
        return [fq1gz] if isSingleEnd else [fq1gz, fq2gz]
    else:
        return [fq1] if isSingleEnd else [fq1, fq2]

if check_sra():
    rule check_sra_strandedness:
        '''Saves the strandedness of the SRAs'''
        output: temp("{dir}/{sra}.spot")
        benchmark: "{dir}/{sra}.spot.benchmark"
        log: "{dir}/{sra}.spot.log"
        shell: "fastq-dump -X 1 -Z --split-spot {wildcards.sra} | wc -l > {output}"

    rule download_sras: # in the future, could use this to check SE vs PE: https://www.biostars.org/p/139422/
        '''Download fastqs from GEO SRA for quality control and alignment'''
        input: "{dir}/{sra}.spot"
        output: "{dir}/{sra}.outlist" # independent of pe/se
        params: fqout=lambda wildcards: sra_output(wildcards)
        benchmark: "{dir}/{sra}.benchmark"
        log: "{dir}/{sra}.log"
        threads: 4
        shell:
            "fasterq-dump -b 10MB -c 100MB -m 1000MB -p --threads {threads}" # use 10x the default memory allocation for larger SRAs
            " --split-files --temp {wildcards.dir} --outdir {wildcards.dir} {wildcards.sra} && "
            "ls {params.fqout} > {output} 2> {log}"
else:
    rule expand_fastqs:
        '''Prepare compressed input fastqs for quality control'''
        input: lambda wildcards: fastqs(wildcards, True)
        output: "{dir}/{fq}.outlist" # independent of pe/se
        params: fqout=lambda wildcards: temp(fastqs(wildcards, False))
        shell: "gunzip -k {input} && ls {params.fqout} > {output}"

def fastp_output(wildcards):
    fastq_list=sra_output(wildcards) if check_sra() else fastqs(wildcards, False)
    fqs_out = [f"{prefix}.trim_{i+1}.fastq.gz" for i in range(len(fastq_list))]
    return fqs_out

rule fastp:
    '''Trim adapters, read quality filtering, make QC outputs'''
    input: "{dir}/{sra}.outlist" if check_sra() else "{dir}/{fq}.outlist"
    output:
        fqlist="{dir}/{sra}.trim.outlist" if check_sra() else "{dir}/{fq}.trim.outlist", # independent of pe/se
        html="{dir}/{sra}.trim.html" if check_sra() else "{dir}/{fq}.trim.html",
        json="{dir}/{sra}.trim.json" if check_sra() else "{dir}/{fq}.trim.json",
    threads: 6
    log: "{dir}/{sra}.trim.log" if check_sra() else "{dir}/{fq}.trim.log"
    params:
        quality=20,
        title="{sra}" if check_sra() else "{fq}",
        outfq=lambda wildcards: fastp_output(wildcards),
        infq=lambda wildcards: fastp_input(wildcards)
    shell:
        "fastp -q {params.quality} " + \
        "-i {params.infq[0]} -o {params.outfq[0]} " if len("{params.infq}".split()) == 1 else \
        "-i {params.infq[0]} -I {params.infq[1]} -o {params.outfq[0]} -O {params.outfq[1]} " + \
        "-h {output.html} -j {output.json} " + \
        "-w {threads} -R {params.title} --detect_adapter_for_pe && " + \
        "ls {params.fqout} > {output.fqlist} &> {log}"

rule hisat2_align_bam:
    '''Align trimmed reads'''
    input:
        f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.1.ht2",
        fq="{dir}/{sra}.trim.outlist" if check_sra() else "{dir}/{fq}.trim.outlist",
        ss=f"data/ensembl/{REF}.{config['release']}.splicesites.txt"
    output:
        "{dir}/align/{sra}.sorted.bam" if check_sra() else "{dir}/align/{fq}.sorted.bam",
    threads: 12
    params:
        compression="9",
        tempprefix="{dir}/align/{sra}.sorted" if check_sra() else "{dir}/align/{fq}.sorted",
        infq=lambda wildcards: fastp_output(wildcards)
    log: "{dir}/align/{sra}.hisat2.log" if check_sra() else "{dir}/align/{fq}.hisat2.log"
    shell:
        "(hisat2 -p {threads} -x data/ensembl/" + REF + ".dna.primary_assembly.karyotypic " + \
        "-U {params.infq[0]} " if len("{params.infq}".split()) == 1 else "-1 {params.infq[0]} -2 {params.infq[1]} " + \
        "--known-splicesite-infile {input.ss} | " # align the suckers
        "samtools view -h -F4 - | " # get mapped reads only
        "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} -) 2> {log} && " # sort them
        "samtools index {output}"

rule hisat2_merge_bams:
    '''Merge the BAM files for each sample'''
    input:
        bams=expand("{{dir}}/align/{sra}.sorted.bam", sra=config["sra"]) if check_sra() else expand("{{dir}}/align/{fq}.sorted.bam", fq=config["fq"])
    output:
        sorted="{dir}/align/combined.sorted.bam",
        stats="{dir}/align/combined.sorted.stats"
    params:
        compression="9",
        tempprefix="{dir}/align/combined.sorted"
    log: "{dir}/align/combined.sorted.log"
    threads: 12
    resources: mem_mb=16000
    shell:
        "(ls {input.bams} | "
        "{{ read firstbam; "
        "samtools view -h ""$firstbam""; "
        "while read bam; do samtools view ""$bam""; done; }} | "
        "samtools view -ubS - | "
        "samtools sort -@ {threads} -l {params.compression} -T {params.tempprefix} -o {output.sorted} - && "
        "samtools index {output.sorted} && "
        "samtools flagstat -@ {threads} {output.sorted} > {output.stats}) 2> {log}"
