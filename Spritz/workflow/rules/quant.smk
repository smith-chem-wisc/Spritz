if check('sra'):
    rule quantify_transcripts_sra:
        input:
            bam="{dir}/align/{sra}.sorted.bam",
            gff="{dir}/isoforms/combined.gtf",
        output:
            gtf=temp("{dir}/isoforms/{sra}.sra.sorted.quant.gtf"),
            gtfgz="{dir}/isoforms/{sra}.sra.sorted.quant.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{sra}.sra.sorted.quant.gtf.log"
        benchmark: "{dir}/isoforms/{sra}.sra.sorted.quant.gtf.benchmark"
        conda: "environments/isoforms_stringtie.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -eB && "
            "gzip -k {output.gtf} 2> {log}"

    rule quantify_ref_transcripts_sra:
        input:
            bam="{dir}/align/{sra}.sorted.bam",
            gff=GFF3,
        output:
            gtf=temp("{dir}/isoforms/{sra}.sra.sorted.quant_ref.gtf"),
            gtfgz="{dir}/isoforms/{sra}.sra.sorted.quant_ref.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{sra}.sra.sorted.quant_ref.gtf.log"
        benchmark: "{dir}/isoforms/{sra}.sra.sorted.quant_ref.gtf.benchmark"
        conda: "environments/isoforms_stringtie.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -eB && "
            "gzip -k {output.gtf} 2> {log}"

if check('sra_se'):
    rule quantify_transcripts_sra_se:
        input:
            bam="{dir}/align/{sra_se}.sorted.bam",
            gff="{dir}/isoforms/combined.gtf",
        output:
            gtf=temp("{dir}/isoforms/{sra_se}.sra_se.sorted.quant.gtf"),
            gtfgz="{dir}/isoforms/{sra_se}.sra_se.sorted.quant.gtf.gz"
        threads: 4
        log: "data/isoforms/{sra_se}.sra_se.sorted.quant.gtf.log"
        benchmark: "data/isoforms/{sra_se}.sra_se.sorted.quant.gtf.benchmark"
        conda: "environments/isoforms_stringtie.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -eB && "
            "gzip -k {output.gtf} 2> {log}"

    rule quantify_ref_transcripts_sra_se:
        input:
            bam="{dir}/align/{sra_se}.sorted.bam",
            gff=GFF3,
        output:
            gtf=temp("{dir}/isoforms/{sra_se}.sra_se.sorted.quant_ref.gtf"),
            gtfgz="{dir}/isoforms/{sra_se}.sra_se.sorted.quant_ref.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{sra_se}.sra_se.sorted.quant_ref.gtf.log"
        benchmark: "{dir}/isoforms/{sra_se}.sra_se.sorted.quant_ref.gtf.benchmark"
        conda: "environments/isoforms_stringtie.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -eB && "
            "gzip -k {output.gtf} 2> {log}"

if check('fq'):
    rule quantify_transcripts_fq:
        input:
            bam="{dir}/align/{fq}.sorted.bam",
            gff="{dir}/isoforms/combined.gtf",
        output:
            gtf=temp("{dir}/isoforms/{fq}.fq.sorted.quant.gtf"),
            gtfgz="{dir}/isoforms/{fq}.fq.sorted.quant.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{fq}.fq.sorted.quant.gtf.log"
        benchmark: "{dir}/isoforms/{fq}.fq.sorted.quant.gtf.benchmark"
        conda: "environments/isoforms_stringtie.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -eB && "
            "gzip -k {output.gtf} 2> {log}"

    rule quantify_ref_transcripts_fq:
        input:
            bam="{dir}/align/{fq}.sorted.bam",
            gff=GFF3,
        output:
            gtf=temp("{dir}/isoforms/{fq}.fq.sorted.quant_ref.gtf"),
            gtfgz="{dir}/isoforms/{fq}.fq.sorted.quant_ref.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{fq}.fq.sorted.quant_ref.gtf.log"
        benchmark: "{dir}/isoforms/{fq}.fq.sorted.quant_ref.gtf.benchmark"
        conda: "environments/isoforms_stringtie.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -eB && "
            "gzip -k {output.gtf} 2> {log}"

if check('fq_se'):
    rule quantify_transcripts_fq_se:
        input:
            bam="{dir}/align/{fq_se}.sorted.bam",
            gff="{dir}/isoforms/combined.gtf",
        output:
            gtf=temp("{dir}/isoforms/{fq_se}.fq_se.sorted.quant.gtf"),
            gtfgz="{dir}/isoforms/{fq_se}.fq_se.sorted.quant.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{fq_se}.fq_se.sorted.quant.gtf.log"
        benchmark: "{dir}/isoforms/{fq_se}.fq_se.sorted.quant.gtf.benchmark"
        conda: "environments/isoforms_stringtie.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -eB && "
            "gzip -k {output.gtf} 2> {log}"

    rule quantify_ref_transcripts_fq_se:
        input:
            bam="{dir}/align/{fq_se}.sorted.bam",
            gff=GFF3,
        output:
            gtf=temp("{dir}/isoforms/{fq_se}.fq_se.sorted.quant_ref.gtf"),
            gtfgz="{dir}/isoforms/{fq_se}.fq_se.sorted.quant_ref.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{fq_se}.fq_se.sorted.quant_ref.gtf.log"
        benchmark: "{dir}/isoforms/{fq_se}.fq_se.sorted.quant_ref.gtf.benchmark"
        conda: "environments/isoforms_stringtie.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -eB && "
            "gzip -k {output.gtf} 2> {log}"

rule make_isoform_quant_dataframe_custom:
    '''Take the stringtie quant results and put them into a combined dataframe'''
    input:
        lambda w:
            ([] if not check('sra') else expand("{{dir}}/isoforms/{sra}.sra.sorted.quant.gtf", sra=config["sra"])) + \
            ([] if not check('sra_se') else expand("{{dir}}/isoforms/{sra_se}.sra_se.sorted.quant.gtf", sra_se=config["sra_se"])) + \
            ([] if not check('fq') else expand("{{dir}}/isoforms/{fq}.fq.sorted.quant.gtf", fq=config["fq"])) + \
            ([] if not check('fq_se') else expand("{{dir}}/isoforms/{fq_se}.fq_se.sorted.quant.gtf", fq_se=config["fq_se"]))
    output: "{dir}/final/isoform_quant.tpms.csv"
    log: "{dir}/final/isoform_quant.tpms.log"
    benchmark: "{dir}/final/isoform_quant.tpms.benchmark"
    conda: "environments/basic.yaml"
    shell:
        "python scripts/SummarizeQuant.py {output} {input}"

rule make_isoform_quant_dataframe_ref:
    '''Take the stringtie quant results and put them into a combined dataframe'''
    input:
        lambda w:
            ([] if not check('sra') else expand("{{dir}}/isoforms/{sra}.sra.sorted.quant_ref.gtf", sra=config["sra"])) + \
            ([] if not check('sra_se') else expand("{{dir}}/isoforms/{sra_se}.sra_se.sorted.quant_ref.gtf", sra_se=config["sra_se"])) + \
            ([] if not check('fq') else expand("{{dir}}/isoforms/{fq}.fq.sorted.quant_ref.gtf", fq=config["fq"])) + \
            ([] if not check('fq_se') else expand("{{dir}}/isoforms/{fq_se}.fq_se.sorted.quant_ref.gtf", fq_se=config["fq_se"]))
    output: "{dir}/final/reference_quant.tpms.csv"
    log: "{dir}/final/reference_quant.tpms.log"
    benchmark: "{dir}/final/reference_quant.tpms.benchmark"
    conda: "environments/basic.yaml"
    shell:
        "python scripts/SummarizeQuant.py {output} {input}"
