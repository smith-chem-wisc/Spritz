## For generating smaller human genomes
rule filter_gff3:
    input: "data/ensembl/" + REF + "." + config["release"] + ".gff3"
    output: "data/ensembl/202122.gff3"
    shell: "grep \"^#\|20\|^21\|^22\" \"data/ensembl/" + REF + "." + config["release"] + ".gff3\" > \"data/ensembl/202122.gff3\""

rule fix_gff3_for_rsem:
    '''This script changes descriptive notes in column 4 to "gene" if a gene row, and it also adds ERCCs to the gene model'''
    input: "data/ensembl/" + REF + "." + config["release"] + ".gff3"
    output: "data/ensembl/" + REF + "." + config["release"] + ".gff3" + ".fix.gff3"
    shell: "python scripts/fix_gff3_for_rsem.py {input} {output}"

rule filter_fa:
    input: "data/ensembl/" + REF + ".dna.primary_assembly.fa"
    output: "data/ensembl/202122.fa"
    script: "../scripts/filter_fasta.py"

rule simulate_variants:
    input: FA
    output: "data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.vcf"
    benchmark: "data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.vcf.benchmark"
    log: "data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.vcf.log"
    shell: "mason_variator -ir {input} -ov {output} 2> {log}"

rule generate_fastqs:
    input:
        fa=FA,
        vcf="data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.vcf"
    output:
        fq1="data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test_1.fastq",
        fq2="data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test_2.fastq",
    benchmark: "data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.benchmark"
    log: "data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.test.log"
    shell: "mason_simulator -ir {input.fa} -n 100000 -iv {input.vcf} -o {output.fq1} -or {output.fq2} 2> {log}"
