REF=config["species"] + "." + config["genome"]

rule download_ensembl_references:
    output:
        "data/ensembl/" + REF + ".dna.primary_assembly.fa",
        "data/ensembl/" + REF + "." + config["release"] + ".gff3",
        "data/ensembl/" + REF + ".pep.all.fa",
        "data/ensembl/" + config["species"] + ".vcf",
    log: "data/ensembl/downloads.log"
    shell:
        "(python scripts/download_ensembl.py " + REF + " && "
        "gunzip data/ensembl/" + REF + ".dna.primary_assembly.fa.gz && "
        "gunzip data/ensembl/" + REF + "." + config["release"] + ".gff3.gz && "
        "gunzip data/ensembl/" + REF + ".pep.all.fa.gz && "
        "gunzip data/ensembl/" + config["species"] + ".vcf.gz) 2> {log}"

rule clean_and_index_vcf:
    input:
        "data/ensembl/" + config["species"] + ".vcf",
    output:
        clean_vcf=temp("data/ensembl/" + config["species"] + ".clean.vcf"),
        vcf="data/ensembl/" + config["species"] + ".ensembl.vcf",
        idx="data/ensembl/" + config["species"] + ".ensembl.vcf.idx"
    log: "data/ensembl/" + config["species"] + ".clean.log"
    shell:
        "(python scripts/clean_vcf.py && "
        "mv {output.clean_vcf} {output.vcf} &&"
        "gatk IndexFeatureFile -F {output.vcf}) 2> {log}"

rule download_chromosome_mappings:
    output: "ChromosomeMappings/" + config["genome"] + "_UCSC2ensembl.txt"
    shell: "rm -rf ChromosomeMappings && git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule reorder_genome_fasta:
    input: "data/ensembl/" + REF + ".dna.primary_assembly.fa"
    output: "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa"
    script: "../scripts/karyotypic_order.py"

rule dict_fa:
    input: "data/ensembl/" + config["species"] + "." + config["genome"] + ".dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/" + config["species"] + "." + config["genome"] + ".dna.primary_assembly.karyotypic.dict"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"

rule tmpdir:
    output:
        temp(directory("tmp")),
        temp(directory("temporary")),
    shell:
        "mkdir tmp && mkdir temporary"

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

rule download_sras:
    output:
        temp("{dir}/{sra,[A-Z0-9]+}_1.fastq"), # constrain wildcards, so it doesn't soak up SRR######.trim_1.fastq
        temp("{dir}/{sra,[A-Z0-9]+}_2.fastq")
    log: "{dir}/{sra}.log"
    threads: 4
    shell:
        "fasterq-dump -p --threads {threads} --split-files --temp {wildcards.dir} --outdir {wildcards.dir} {wildcards.sra} 2> {log}"

rule compress_fastqs:
    input:
        temp("{dir}/{sra,[A-Z0-9]+}_1.fastq"),
        temp("{dir}/{sra,[A-Z0-9]+}_2.fastq")
    output:
        temp("{dir}/{sra,[A-Z0-9]+}_1.fastq.gz"),
        temp("{dir}/{sra,[A-Z0-9]+}_2.fastq.gz")
    threads: 2
    shell:
        "gzip {input[0]} & gzip {input[1]}"
