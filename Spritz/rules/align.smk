import os

REF=config["species"] + "." + config["genome"]

rule directories:
    output: directory("data/ensembl/{REF}.dna.primary_assembly.karyotypic/")
    shell: "mkdir -p data/ensembl/{REF}.dna.primary_assembly.karyotypic/"

rule hisat_genome:
    input:
        fa="data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        gtf="data/ensembl/{REF}." + config["release"] + ".gff3",
    threads: 12
    output:
        idx="data/ensembl/{REF}.dna.primary_assembly.karyotypic.1.ht2",
        finished="data/ensembl/done_building_hisat_genome{REF}.txt",
    benchmark: "data/ensembl/{REF}.hisatbuild.benchmark"
    log: "data/ensembl/{REF}.hisatbuild.log"
    shell: 
        "(hisat2-build -p {threads} data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa"
        " data/ensembl/{REF}.dna.primary_assembly.karyotypic && touch {output.finished}) &> {log}"

rule hisat2_splice_sites:
    input: "data/ensembl/" + REF + "." + config["release"] + ".gff3"
    output: "data/ensembl/" + REF + "." + config["release"] + ".splicesites.txt"
    shell: "hisat2_extract_splice_sites.py {input} > {output}"

def input_fq_args(fastqs):
    fqs=fastqs.split()
    if len(fqs) == 1:
        return f"-U {fqs[0]}"
    else:
        return f"-1 {fqs[0]} -2 {fqs[1]}"

rule fastp:
    '''Trim adapters, read quality filtering, make QC outputs'''
    input:
        fq1="{dir}/{sra}_1.fastq" if check_sra() else "{dir}/{fq}_1.fastq",
        fq2="{dir}/{sra}_2.fastq" if check_sra() else "{dir}/{fq}_2.fastq",
    output:
        fq1="{dir}/{sra}.trim_1.fastq.gz" if check_sra() else "{dir}/{fq}.trim_1.fastq.gz",
        fq2="{dir}/{sra}.trim_2.fastq.gz" if check_sra() else "{dir}{fq}.trim_2.fastq.gz",
        html="{dir}/{sra}.trim.html" if check_sra() else "{dir}/{fq}.trim.html",
        json="{dir}/{sra}.trim.json" if check_sra() else "{dir}/{fq}.trim.json",
    threads: 6
    log: "{dir}/{sra}.trim.log" if check_sra() else "{dir}/{fq}.trim.log"
    params:
        quality=20,
        title="{sra}" if check_sra() else "{fq}"
    shell:
        "fastp -q {params.quality} -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} "
        "-h {output.html} -j {output.json} "
        "-w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

rule hisat2_align_bam:
    '''Align trimmed reads'''
    input:
        "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.1.ht2",
        fq1="{dir}/{sra}.trim_1.fastq.gz" if check_sra() else "{dir}/{fq}.trim_1.fastq.gz",
        fq2="{dir}/{sra}.trim_2.fastq.gz" if check_sra() else "{dir}{fq}.trim_2.fastq.gz",
        ss="data/ensembl/" + REF + "." + config["release"] + ".splicesites.txt"
    output:
        sorted="{dir}/{sra}.sorted.bam" if check_sra() else "{dir}/{fq}.sorted.bam",
    threads: 12
    params:
        compression="9",
        tempprefix="{dir}/{sra}.sorted" if check_sra() else "{dir}/{fq}.sorted",
    log: "{dir}/{sra}.hisat2.log" if check_sra() else "{dir}/{fq}.hisat2.log"
    shell:
        "(hisat2 -p {threads} -x data/ensembl/" + REF + ".dna.primary_assembly.karyotypic -1 {input.fq1} -2 {input.fq2} --known-splicesite-infile {input.ss} | " # align the suckers
        "samtools view -h -F4 - | " # get mapped reads only
        "samtools sort -l {params.compression} -T {params.tempprefix} -o {output.sorted} -) 2> {log} && " # sort them
        "samtools index {output}"

rule hisat2_merge_bams:
    '''Merge the BAM files for each sample'''
    input:
        bams=expand("{{dir}}/{sra}.sorted.bam", sra=config["sra"]) if check_sra() else expand("{{dir}}/{fq}.sorted.bam", fq=config["fq"])
    output:
        sorted="{dir}/combined.sorted.bam",
        stats="{dir}/combined.sorted.stats"
    params:
        compression="9",
        tempprefix="{dir}/combined.sorted"
    log: "{dir}/combined.sorted.log"
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
