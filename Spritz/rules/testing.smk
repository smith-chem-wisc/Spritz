rule filter_gff3:
    '''Generating a smaller GFF3 file'''
    input: f"data/ensembl/{REF}/{config['release']}.gff3"
    output: "data/ensembl/202122.gff3"
    params:
        ref=REF,
        release=config['release']
    log: "data/ensembl/202122.gff3.log"
    shell: "grep \"^#\|20\|^21\|^22\" \"data/ensembl/{params.ref}.{params.release}.gff3\" > \"data/ensembl/202122.gff3\" 2> {log}"

rule filter_fa:
    '''Generating a smaller genome fasta file'''
    input: f"data/ensembl/{REF}.dna.primary_assembly.fa"
    output: "data/ensembl/202122.fa"
    log: "data/ensembl/202122.fa.log"
    script: "../scripts/filter_fasta.py 2> {log}"

rule fix_gff3_for_rsem:
    '''This script changes descriptive notes in column 4 to "gene" if a gene row, and it also adds ERCCs to the gene model'''
    input: f"data/ensembl/{REF}/{config['release']}.gff3"
    output: f"data/ensembl/{REF}/{config['release']}.gff3.fix.gff3"
    log: f"data/ensembl/{REF}/{config['release']}.gff3.fix.log"
    shell: "python scripts/fix_gff3_for_rsem.py {input} {output} 2> {log}"

rule simulate_variants:
    input: FA
    output: f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.vcf"
    benchmark: f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.vcf.benchmark"
    log: f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.vcf.log"
    shell: "mason_variator -ir {input} -ov {output} 2> {log}"

rule generate_fastqs:
    input:
        fa=FA,
        vcf=f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.vcf"
    output:
        fq1=f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test_1.fastq",
        fq2=f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test_2.fastq",
    benchmark: f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.benchmark"
    log: f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.log"
    shell: "mason_simulator -ir {input.fa} -n 100000 -iv {input.vcf} -o {output.fq1} -or {output.fq2} 2> {log}"
