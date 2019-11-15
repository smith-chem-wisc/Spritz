rule download_ensembl_references:
    output:
        gfa=GENOME_FA,
        gff=ENSEMBL_GFF,
        pfa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".pep.all.fa",
        vcf="data/ensembl/common_all_20170710.vcf",
        vcfidx="data/ensembl/common_all_20170710.vcf.idx"
    log: "data/ensembl/downloads.log"
    shell:
        "if [ -e {input.gfa}.gz ] && [ -e {input.gff}.gz ] && [ -e {input.pfa}.gz ] && [ -e {input.vcf}.gz ]; then; "
            "gunzip -c {input.gfa} > {output.gfa} && "
            "gunzip -c {input.gff} > {output.gff} && "
            "gunzip -c {input.pfa} > {output.pfa} && "
            "gunzip -c {input.vcf} > {output.vcf} && "
            "gatk IndexFeatureFile -F {output.vcf}; "
        "else; "
            "(wget -O - ftp://ftp.ensembl.org/pub/release-" + ENSEMBL_VERSION + "//fasta/homo_sapiens/dna/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.fa.gz | "
            "gunzip -c > {output.gfa} && "
            "wget -O - ftp://ftp.ensembl.org/pub/release-" + ENSEMBL_VERSION + "/gff3/homo_sapiens/Homo_sapiens." + GENEMODEL_VERSION + ".gff3.gz | "
            "gunzip -c > {output.gff} && "
            "wget -O - ftp://ftp.ensembl.org/pub/release-" + ENSEMBL_VERSION + "//fasta/homo_sapiens/pep/Homo_sapiens." + GENOME_VERSION + ".pep.all.fa.gz | "
            "gunzip -c > {output.pfa} && "
            "wget -O - ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz | "
            "gunzip -c > {output.vcf} && "
            "gatk IndexFeatureFile -F {output.vcf}) 2> {log}; "
        "fi"

rule download_chromosome_mappings:
    output: "ChromosomeMappings/GRCh38_UCSC2ensembl.txt"
    shell: "git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule reorder_genome_fasta:
    input: "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.fa"
    output: "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa"
    script: "../scripts/karyotypic_order.py"

rule convert_ucsc2ensembl:
    input:
        "data/ensembl/common_all_20170710.vcf",
        "ChromosomeMappings/GRCh38_UCSC2ensembl.txt"
    output:
        "data/ensembl/common_all_20170710.ensembl.vcf",
    script:
        "../scripts/convert_ucsc2ensembl.py"

rule index_ucsc2ensembl:
    input: "data/ensembl/common_all_20170710.ensembl.vcf"
    output: "data/ensembl/common_all_20170710.ensembl.vcf.idx"
    shell: "gatk IndexFeatureFile -F {input}"

rule filter_gff3:
    input: ENSEMBL_GFF
    output: "data/ensembl/202122.gff3"
    shell: "grep \"^#\|^20\|^21\|^22\" \"data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".gff3\" > \"data/ensembl/202122.gff3\""

rule fix_gff3_for_rsem:
    '''This script changes descriptive notes in column 4 to "gene" if a gene row, and it also adds ERCCs to the gene model'''
    input: ENSEMBL_GFF
    output: ENSEMBL_GFF + ".fix.gff3"
    shell: "python scripts/fix_gff3_for_rsem.py {input} {output}"

rule filter_fa:
    input: "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.fa"
    output: "data/ensembl/202122.fa"
    script: "../scripts/filter_fasta.py"

rule download_sras:
    output:
        temp("data/{sra,[A-Z0-9]+}_1.fastq"), # constrain wildcards, so it doesn't soak up SRR######.trim_1.fastq
        temp("data/{sra,[A-Z0-9]+}_2.fastq")
    log: "data/{sra}.log"
    threads: 4
    shell:
        "fasterq-dump --progress --threads {threads} --split-files --outdir data {wildcards.sra} 2> {log}"

rule compress_fastqs:
    input:
        temp("data/{sra,[A-Z0-9]+}_1.fastq"),
        temp("data/{sra,[A-Z0-9]+}_2.fastq")
    output:
        temp("data/{sra,[A-Z0-9]+}_1.fastq.gz"),
        temp("data/{sra,[A-Z0-9]+}_2.fastq.gz")
    threads: 2
    shell:
        "gzip {input[0]} & gzip {input[1]}"
