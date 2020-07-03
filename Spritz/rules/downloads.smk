REF=config["species"] + "." + config["genome"]
SPECIES_LOWER = config["species"].lower()

rule download_ensembl_references:
    output:
        gfa="data/ensembl/" + REF + ".dna.primary_assembly.fa",
        gff3="data/ensembl/" + REF + "." + config["release"] + ".gff3",
        pfa="data/ensembl/" + REF + ".pep.all.fa",
    benchmark: "data/ensembl/downloads.benchmark"
    log: "data/ensembl/downloads.log"
    shell:
        "(python scripts/download_ensembl.py {REF}.{ENSEMBL_VERSION} not && "
        "gunzip {output.gfa}.gz {output.gff3}.gz {output.pfa}.gz) 2> {log}"

if SPECIES_LOWER == "homo_sapiens":
    rule download_dbsnp_vcf:
        '''Download dbsnp known variant sites if we are analyzing human data'''
        input: "ChromosomeMappings/" + config["genome"] + "_UCSC2ensembl.txt"
        output: "data/ensembl/" + config["species"] + ".ensembl.vcf",
        benchmark: "data/ensembl/downloads_dbsnp_vcf.benchmark"
        log: "data/ensembl/downloads_dbsnp_vcf.log"
        shell:
            "(wget -O - https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz | "
            "zcat - | python scripts/convert_ucsc2ensembl.py > {output}) 2> {log}"
else:
    # first get the possible VCF urls; note that Ensembl has started listing variants for each chromosome for human but not other species, but that may change
    vcf1 = f"http://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/variation/vcf/{SPECIES_LOWER}/{SPECIES}.vcf.gz"
    vcf2 = f"http://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/variation/vcf/{SPECIES_LOWER}/{SPECIES_LOWER}.vcf.gz"

    rule download_ensembl_vcf:
        '''Use Ensembl known variant sites if we are analyzing nonhuman data'''
        output: "data/ensembl/" + config["species"] + ".ensembl.vcf",
        benchmark: "data/ensembl/downloads_ensembl_vcf.benchmark"
        log: "data/ensembl/downloads_ensembl_vcf.log"
        shell: "((wget -O - {vcf1} || wget -O - {vcf2}) | zcat - | python scripts/clean_vcf.py > {output}) 2> {log}"

rule index_ensembl_vcf:
    input: "data/ensembl/" + config["species"] + ".ensembl.vcf"
    output: "data/ensembl/" + config["species"] + ".ensembl.vcf.idx"
    log: "data/ensembl/" + config["species"] + ".ensembl.vcf.idx.log"
    shell: "gatk IndexFeatureFile -F {input} 2> {log}"

rule download_chromosome_mappings:
    output: "ChromosomeMappings/" + config["genome"] + "_UCSC2ensembl.txt"
    log: "ChromosomeMappings/download_chromosome_mappings.log"
    shell:
        "(if [ -d ChromosomeMappings ]; then rm -rf ChromosomeMappings; fi && "
        "git clone https://github.com/dpryan79/ChromosomeMappings.git) 2> {log}"

rule reorder_genome_fasta:
    input: "data/ensembl/" + REF + ".dna.primary_assembly.fa"
    output: "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa"
    benchmark: "data/ensembl/karyotypic_order.benchmark"
    log: "data/ensembl/karyotypic_order.log"
    shell: "python scripts/karyotypic_order.py 2> {log}"

rule dict_fa:
    input: "data/ensembl/" + config["species"] + "." + config["genome"] + ".dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/" + config["species"] + "." + config["genome"] + ".dna.primary_assembly.karyotypic.dict"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"

rule download_sras:
    output:
        temp("{dir}/{sra,[A-Z0-9]+}_1.fastq"), # constrain wildcards, so it doesn't soak up SRR######.trim_1.fastq
        temp("{dir}/{sra,[A-Z0-9]+}_2.fastq")
    benchmark: "{dir}/{sra}.benchmark"
    log: "{dir}/{sra}.log"
    threads: 4
    shell:
        "fasterq-dump -p --threads {threads} --split-files --temp {wildcards.dir} --outdir {wildcards.dir} {wildcards.sra} 2> {log}"

rule tmpdir:
    output:
        temp(directory("tmp")),
        temp(directory("temporary")),
    shell:
        "mkdir tmp && mkdir temporary"
