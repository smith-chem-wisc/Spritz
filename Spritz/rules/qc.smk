rule download_adapters:
    output:
        temp(directory("BBMap")),
        "data/qc/adapters.fa"
    shell: "git clone --depth 1 https://github.com/BioInfoTools/BBMap.git && cp BBMap/resources/adapters.fa data/qc"

rule skewer:
    input:
        fq1="data/{sra}_1.fastq.gz",
        fq2="data/{sra}_2.fastq.gz",
        adapters="data/qc/adapters.fa"
    output:
        fq1="data/trimmed{q}/{sra}.trim_1.fastq.gz",
        fq2="data/trimmed{q}/{sra}.trim_2.fastq.gz",
    threads: 12
    log: "data/trimmed{q}/{sra}-trimmed.status"
    params:
        quality="{q}",
        ext="data/trimmed{q}/{sra}"
    shell:
        "skewer -q {params.quality} -o {params.ext}"
        " -t {threads} -x {input.adapters} {input.fq1} {input.fq2} &> {log} && "
        "mv {params.ext}-trimmed-pair1.fastq {params.ext}.trim_1.fastq && "
        "mv {params.ext}-trimmed-pair2.fastq {params.ext}.trim_2.fastq && "
        "gzip {params.ext}.trim_*.fastq"

rule fastqc_analysis:
    input:
        fq1=["data/{sra}_1.fastq.gz", "data/trimmed{q}/{sra}.trim_1.fastq.gz"],
        fq2=["data/{sra}_2.fastq.gz", "data/trimmed{q}/{sra}.trim_2.fastq.gz"],
    output:
        fq1=["data/trimmed{q}/{sra}_1_fastqc.html", "data/trimmed{q}/{sra}_1_fastqc.zip", "data/trimmed{q}/{sra}.trim_1_fastqc.html", "data/trimmed{q}/{sra}.trim_1_fastqc.zip"],
        fq2=["data/trimmed{q}/{sra}_2_fastqc.html", "data/trimmed{q}/{sra}_2_fastqc.zip", "data/trimmed{q}/{sra}.trim_2_fastqc.html", "data/trimmed{q}/{sra}.trim_2_fastqc.zip"],
    log: "data/trimmed{q}/{sra}.fastqc.log"
    threads: 6
    shell:
        "fastqc -o data/trimmed{q} -t {threads} {input.fq1} {input.fq2} 2> {log}"
