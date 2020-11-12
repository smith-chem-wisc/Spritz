rule rsem_star_genome:
    '''Create an RSEM reference with STAR indices'''
    input:
        gfa=f"data/ensembl/{REF}.dna.primary_assembly.fa",
        gff=f"data/ensembl/{REF}.{config['release']}.gff3.fix.gff3"
    output:
        f"{REFSTAR_PREFIX}.gtf",
        f"{REFSTAR_FOLDER}SA",
    params: refstar_prefix=lambda w, output: os.path.splitext(output[0])[0]
    threads: 99
    resources: mem_mb=60000
    log: "data/ensembl/prepare-reference.log"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star "
        "--gff3 {input.gff} \"{input.gfa}\" {params.refstar_prefix}) 2> {log}"

if check('sra'):
    rule rsem_star_align_sra:
        '''Align to transcripts with STAR and quantify with RSEM'''
        input:
            gtf=f"{REFSTAR_PREFIX}.gtf",
            suffix=f"{REFSTAR_FOLDER}SA",
            fq1="{dir}/trimmed/{sra}.trim_1.fastq.gz",
            fq2="{dir}/trimmed/{sra}.trim_2.fastq.gz",
        output:
            "{dir}/quant/{sra}.sra.isoforms.results",
            "{dir}/quant/{sra}.sra.genes.results",
            "{dir}/quant/{sra}.sra.time",
            directory("{dir}/quant/{sra}.sra.stat"),
        params: refstar_prefix=lambda w, input: os.path.splitext(input[0])[0]
        resources: mem_mb=50000
        threads: 12
        log: "{dir}/{sra}.sra.calculate-expression.log"
        shell:
            "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" not doing credibility intervals for now; they take a long time to calc.
            " --num-threads {threads} --paired-end <(zcat {input.fq1}) <(zcat {input.fq2}) "
            "{params.refstar_prefix} {wildcards.dir}/quant/{wildcards.sra}) &> {log}"
if check('sra_se'):
    rule rsem_star_align_sra_se:
        '''Align to transcripts with STAR and quantify with RSEM'''
        input:
            gtf=f"{REFSTAR_PREFIX}.gtf",
            suffix=f"{REFSTAR_FOLDER}SA",
            fq1="{dir}/trimmed/{sra_se}.trim.fastq.gz",
        output:
            "{dir}/quant/{sra_se}.sra_se.isoforms.results",
            "{dir}/quant/{sra_se}.sra_se.genes.results",
            "{dir}/quant/{sra_se}.sra_se.time",
            directory("{dir}/quant/{sra_se}.sra_se.stat"),
        params: refstar_prefix=lambda w, input: os.path.splitext(input[0])[0]
        resources: mem_mb=50000
        threads: 12
        log: "{dir}/{sra_se}.sra_se.calculate-expression.log"
        shell:
            "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" not doing credibility intervals for now; they take a long time to calc.
            " --num-threads {threads} <(zcat {input.fq1}) "
            "{params.refstar_prefix} {wildcards.dir}/quant/{wildcards.sra_se}) &> {log}"
if check('fq'):
    rule rsem_star_align_fq:
        '''Align to transcripts with STAR and quantify with RSEM'''
        input:
            gtf=f"{REFSTAR_PREFIX}.gtf",
            suffix=f"{REFSTAR_FOLDER}SA",
            fq1="{dir}/{fq}.trim_1.fastq.gz",
            fq2="{dir}/{fq}.trim_2.fastq.gz",
        output:
            "{dir}/quant/{fq}.fq.isoforms.results",
            "{dir}/quant/{fq}.fq.genes.results",
            "{dir}/quant/{fq}.fq.time",
            directory("{dir}/quant/{fq}.fq.stat"),
        params: refstar_prefix=lambda w, input: os.path.splitext(input[0])[0]
        resources: mem_mb=50000
        threads: 12
        log: "{dir}/{fq}.fq.calculate-expression.log"
        shell:
            "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" not doing credibility intervals for now; they take a long time to calc.
            " --num-threads {threads} --paired-end <(zcat {input.fq1}) <(zcat {input.fq2}) "
            "{params.refstar_prefix} {wildcards.dir}/quant/{wildcards.fq}) &> {log}"

if check('fq_se'):
    rule rsem_star_align_fq_se:
        '''Align to transcripts with STAR and quantify with RSEM'''
        input:
            gtf=f"{REFSTAR_PREFIX}.gtf",
            suffix=f"{REFSTAR_FOLDER}SA",
            fq1="{dir}/{fq_se}.trim_1.fastq.gz",
        output:
            "{dir}/quant/{fq_se}.fq_se.isoforms.results",
            "{dir}/quant/{fq_se}.fq_se.genes.results",
            "{dir}/quant/{fq_se}.fq_se.time",
            directory("{dir}/quant/{fq_se}.fq_se.stat"),
        params: refstar_prefix=lambda w, input: os.path.splitext(input[0])[0]
        resources: mem_mb=50000
        threads: 12
        log: "{dir}/{fq_se}.fq_se.calculate-expression.log"
        shell:
            "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" not doing credibility intervals for now; they take a long time to calc.
            " --num-threads {threads} <(zcat {input.fq1}) "
            "{params.refstar_prefix} {wildcards.dir}/quant/{wildcards.fq_se}) &> {log}"

rule make_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        lambda w:
            [] if not check('sra') else expand("{{dir}}/quant/{sra}.sra.gene.results", sra=config["sra"]) + \
            [] if not check('sra_se') else expand("{{dir}}/quant/{sra_se}.sra_se.gene.results", sra_se=config["sra_se"]) + \
            [] if not check('fq') else expand("{{dir}}/quant/{fq}.fq.gene.results", fq=config["fq"]) + \
            [] if not check('fq_se') else expand("{{dir}}/quant/{fq_se}.fq_se.gene.results", fq_se=config["fq_se"]),
        gff=f"data/ensembl/{REF}.{config['release']}.gff3.fix.gff3"
    output:
        counts="{dir}/Counts.csv",
        names="{dir}/IdsToNames.csv",
        tpms="{dir}/Tpms.csv"
    log: "{dir}/make_rsem_dataframe.log"
    shell:
        "python scripts/make_rsem_dataframe.py {input.gff} {output.counts} {output.tpms} {output.names} 2> {log}"
