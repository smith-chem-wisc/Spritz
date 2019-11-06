rule rsem_star_genome:
    '''Create an RSEM reference with STAR indices'''
    input:
        gfa=FA,
        gff=GFF3 + ".fix.gff3"
    output:
        REFSTAR_PREFIX + ".gtf",
        suffix = REFSTAR_FOLDER + "SA"
    threads: 99
    resources: mem_mb=60000
    log: "data/ensembl/prepare-reference.log"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star --gff3 {input.gff} \"{input.gfa}\" " + REFSTAR_PREFIX +
        ") 2> {log}"

rule rsem_star_align:
    '''Align to transcripts with STAR and quantify with RSEM'''
    input:
        suffix=REFSTAR_FOLDER + "SA",
        gtf=REFSTAR_PREFIX + ".gtf",
        fq1="{dir}/trimmed/{sra}.trim_1.fastq.gz" if check_sra() else "{dir}/{fq}_1.fastq.gz",
        fq2="{dir}/trimmed/{sra}.trim_2.fastq.gz" if check_sra() else "{dir}/{fq}_2.fastq.gz",
    output:
        "{dir}/{sra}.isoforms.results",
        "{dir}/{sra}.genes.results",
        "{dir}/{sra}.time",
        directory("{dir}/{sra}.stat"),
    resources: mem_mb=50000
    threads: 12
    log: "{dir}/{sra}calculate-expression.log"
    shell:
        "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" not doing credibility intervals for now; they take a long time to calc.
        " --num-threads {threads} --paired-end <(zcat {input.fq1}) <(zcat {input.fq2}) " + REFSTAR_PREFIX + " {dir}/{wildcards.sra}) &> {log}"

rule make_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("{dir}/{sra}.genes.results", sra=config["sra"]),
        gff="data/ensembl/" + REF + "." + config["release"] + ".gff3" + ".fix.gff3"
    output:
        counts="{dir}/Counts.csv",
        names="{dir}/IdsToNames.csv",
        tpms="{dir}/Tpms.csv"
    shell:
        "python scripts/make_rsem_dataframe.py {input.gff} {output.counts} {output.tpms} {output.names}"
