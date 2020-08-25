import os

rule directories:
    output: directory("data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic/")
    shell: "mkdir -p data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic/"

rule star_setup:
    output: "STAR-2.6.0c/bin/Linux_x86_64/STAR"
    shell:
        "wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz"
        "tar xvf 2.6.0c.tar.gz"
        "rm 2.6.0c.tar.gz"

rule star_genome_generate:
    input:
        star="STAR-2.7.0e/bin/Linux_x86_64/STAR",
        genomeDir=directory("data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic"),
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        gff="data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".gff3"
    output:
        "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic/SA"
    shell:
        "{input.star} --runMode genomeGenerate --runThreadN {threads} --genomeDir {input.genomeDir} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"

rule hisat_genome:
    input:
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        gtf="data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".gff3"
    threads: 12
    output: "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.1.ht2"
    shell: "hisat2-build -p {threads} data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic"

rule hisat2_splice_sites:
    input: "data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".gff3"
    output: "data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".splicesites.txt"
    shell: "hisat2_extract_splice_sites.py {input} > {output}"

rule hisat2_align_bam:
    input:
        "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.1.ht2",
        # "data/trimmed{q}/{sra}.trim_1_fastqc.html", # trigger QC analysis after trimming
        fq1="data/trimmed{q}/{sra}.trim_1.fastq.gz",
        fq2="data/trimmed{q}/{sra}.trim_2.fastq.gz",
        ss="data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".splicesites.txt"
    output:
        sorted="data/{q}/{sra}.sorted.bam"
    threads: 12
    params:
        compression="9",
        tempprefix="data/{q}/{sra}.sorted"
    log: "data/{q}/{sra}.hisat2.log"
    shell:
        "(hisat2 -p {threads} -x data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic -1 {input.fq1} -2 {input.fq2} --known-splicesite-infile {input.ss} | " # align the suckers
        "samtools view -h -F4 - | " # get mapped reads only
        "samtools sort -l {params.compression} -T {params.tempprefix} -o {output.sorted} -) 2> {log} && " # sort them
        "samtools index {output}"

rule hisat2_merge_bams:
    input:
        bams=expand("data/{{q}}/{sra}.sorted.bam", sra=config["sra"])
    output:
        sorted="data/{q}/combined.sorted.bam",
        stats="data/{q}/combined.sorted.stats"
    params:
        compression="9",
        tempprefix="data/{q}/combined.sorted"
    log: "data/{q}/combined.sorted.log"
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
