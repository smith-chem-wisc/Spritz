import os

REF=config["species"] + "." + config["genome"]

rule directories:
    output: directory("data/ensembl/" + REF + ".dna.primary_assembly.karyotypic/")
    shell: "mkdir -p data/ensembl/" + REF + ".dna.primary_assembly.karyotypic/"

rule star_setup:
    output: "STAR-2.6.0c/bin/Linux_x86_64/STAR"
    shell:
        "wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz"
        "tar xvf 2.6.0c.tar.gz"
        "rm 2.6.0c.tar.gz"

rule star_genome_generate:
    input:
        star="STAR-2.7.0e/bin/Linux_x86_64/STAR",
        genomeDir=directory("data/ensembl/" + REF + ".dna.primary_assembly.karyotypic"),
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        gff="data/ensembl/" + REF + "." + config["release"] + ".gff3"
    output:
        "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic/SA"
    shell:
        "{input.star} --runMode genomeGenerate --runThreadN {threads} --genomeDir {input.genomeDir} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"

rule hisat_genome:
    input:
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        gtf="data/ensembl/" + REF + "." + config["release"] + ".gff3"
    threads: 12
    output: "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.1.ht2"
    shell: "hisat2-build -p {threads} data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa data/ensembl/" + REF + ".dna.primary_assembly.karyotypic"

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

def check_sra():
    docheck = 'sra' in config and config["sra"] is not None) and len(config["sra"]) > 0
    return docheck

rule hisat2_align_bam:
    input:
        "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.1.ht2",
        "{dir}/trimmed/{sra}.trim_1_fastqc.html" if check_sra() else "{dir}/trimmed/{fq}.trim_1_fastqc.html",
        fq1="{dir}/trimmed/{sra}.trim_1.fastq.gz" if check_sra() else "{dir}/trimmed/{fq}.trim_1.fastq.gz",
        fq2="{dir}/trimmed/{sra}.trim_2.fastq.gz" if check_sra() else "{dir}/trimmed/{fq}.trim_2.fastq.gz",
        ss="data/ensembl/" + REF + "." + config["release"] + ".splicesites.txt"
    output:
        sorted="{dir}/{sra}.sorted.bam" if check_sra() else "{dir}/{fq,[A-Z0-9]+}.sorted.bam",
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
