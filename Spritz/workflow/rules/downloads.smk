SPECIES_LOWER = config['species'].lower()
PROTOCOL = "http"

if DIVISION == "bacteria":
    rule download_ensembl_bacteria_metadata:
        '''The species index, which is the only place a species\' collection directory is recorded'''
        output: f"../resources/ensembl/species_EnsemblBacteria.txt"
        params: url=f"https://ftp.ebi.ac.uk/ensemblgenomes/pub/bacteria/release-{ENSEMBL_VERSION}/species_EnsemblBacteria.txt"
        benchmark: "../resources/ensembl/downloads_bacteria_metadata.benchmark"
        log: "../resources/ensembl/downloads_bacteria_metadata.log"
        conda: "../envs/downloads.yaml"
        shell: "wget -O {output} {params.url} 2> {log}"

    rule download_ensembl_bacteria_references:
        '''Genome, gene model and proteome for one bacterial strain from Ensembl Genomes'''
        input: metadata="../resources/ensembl/species_EnsemblBacteria.txt"
        output:
            gfa=GENOME_FA,
            gff3=ENSEMBL_GFF,
            pfa=f"../resources/ensembl/{REF}.pep.all.fa",
        params: species=SPECIES, assembly=GENOME_VERSION, release=ENSEMBL_VERSION
        benchmark: "../resources/ensembl/downloads.benchmark"
        log: "../resources/ensembl/downloads.log"
        conda: "../envs/downloads.yaml"
        # A script rather than wget: the collection directory has to be looked up per species, and
        # the DNA filename carries an unpredictable trailing underscore, so it is read out of the
        # directory's CHECKSUMS index instead of being templated. See scripts/ensembl_bacteria.py.
        shell:
            "python scripts/ensembl_bacteria.py"
            " --species {params.species} --assembly {params.assembly} --release {params.release}"
            " --metadata {input.metadata}"
            " --genome-out {output.gfa} --gff3-out {output.gff3} --protein-out {output.pfa} 2> {log}"

else:
    rule download_ensembl_references:
        output:
            gfa=GENOME_FA,
            gff3=ENSEMBL_GFF,
            pfa=f"../resources/ensembl/{REF}.pep.all.fa",
        params:
            primary=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/dna/{REF}.dna.primary_assembly.fa.gz",
            toplevel=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/dna/{REF}.dna.toplevel.fa.gz",
            gff=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/gff3/{SPECIES_LOWER}/{REF}.{ENSEMBL_VERSION}.gff3.gz",
            gffdir=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/gff3/{SPECIES_LOWER}/",
            gffpattern=f"{REF}\\.[0-9]+\\.gff3\\.gz",
            pep=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/pep/{REF}.pep.all.fa.gz",
        benchmark: "../resources/ensembl/downloads.benchmark"
        log: "../resources/ensembl/downloads.log"
        conda: "../envs/downloads.yaml"
        shell:
            # The gff3 is not always numbered with the Ensembl release: species Ensembl imports from
            # elsewhere carry their own gene set version, so release 116 ships
            # Caenorhabditis_elegans.WBcel235.63.gff3.gz. Read the name off the directory listing and
            # fall back to the assumed one if the listing cannot be read.
            "((wget -O - {params.primary} || wget -O - {params.toplevel}) | gunzip -c - > {output.gfa} && "
            "gff3name=$(wget -qO- {params.gffdir} | grep -oE '{params.gffpattern}' | sort -u | head -1) && "
            "if [ -n \"$gff3name\" ]; then gff3url={params.gffdir}$gff3name; else gff3url={params.gff}; fi && "
            "wget -O - \"$gff3url\" | gunzip -c - > {output.gff3} && "
            "wget -O - {params.pep} | gunzip -c - > {output.pfa}) 2> {log}"

if SPECIES_LOWER == "homo_sapiens":
    rule download_dbsnp_vcf:
        '''Download dbsnp known variant sites if we are analyzing human data'''
        input: f"../resources/ChromosomeMappings/{GENOME_VERSION}_UCSC2ensembl.txt"
        output: f"../resources/ensembl/{SPECIES}.ensembl.vcf",
        params:
            vcf="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz"
        benchmark: "../resources/ensembl/downloads_dbsnp_vcf.benchmark"
        log: "../resources/ensembl/downloads_dbsnp_vcf.log"
        conda: "../envs/downloads.yaml"
        shell:
            "(wget -O - {params.vcf} | zcat - | python scripts/convert_ucsc2ensembl.py > {output}) 2> {log}"

    rule reorder_genome_fasta:
        '''Reorder the ensembl genome to match the dbsnp VCF

        Runs whether or not dbSNP is downloaded: its real job is to produce the karyotypic FASTA that
        the SnpEff database build, samtools faidx and the modification transfer all read.
        '''
        input: GENOME_FA
        output: KARYOTYPIC_GENOME_FA
        benchmark: "../resources/ensembl/karyotypic_order.benchmark"
        log: "../resources/ensembl/karyotypic_order.log"
        conda: "../envs/downloads.yaml"
        shell: "python scripts/karyotypic_order.py 2> {log}"
        
elif DIVISION == "bacteria":
    # Ensembl Bacteria publishes no variant sites - its variation/ directory holds only
    # indexed_vep_cache, with no vcf/ - and GATK base recalibration is the only thing that reads
    # them, so there is no rule to define. This is why a bacterial reference requires a supplied
    # VCF: variants cannot be called from reads without known sites. SpritzCMD rejects the
    # combination up front rather than letting it fail here as a missing input.
    rule rename_genome_fasta_bacteria:
        '''Rename the genome fasta so the other rules work with bacterial genomes'''
        input: GENOME_FA
        output: KARYOTYPIC_GENOME_FA
        benchmark: "../resources/ensembl/karyotypic_order_rename.benchmark"
        log: "../resources/ensembl/karyotypic_order_rename.log"
        conda: "../envs/downloads.yaml"
        shell: "cp {input} {output} 2> {log}"

else:
    rule download_ensembl_vcf:
        '''
        Use Ensembl known variant sites if we are analyzing nonhuman data.
        Note that Ensembl has started listing variants for each chromosome for human, but not other species, but that may change
        '''
        output: f"../resources/ensembl/{SPECIES}.ensembl.vcf",
        params:
            vcf1 = f"http://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/variation/vcf/{SPECIES_LOWER}/{SPECIES}.vcf.gz",
            vcf2 = f"http://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/variation/vcf/{SPECIES_LOWER}/{SPECIES_LOWER}.vcf.gz",
        benchmark: "../resources/ensembl/downloads_ensembl_vcf.benchmark"
        log: "../resources/ensembl/downloads_ensembl_vcf.log"
        conda: "../envs/downloads.yaml"
        shell: "((wget -O - {params.vcf1} || wget -O - {params.vcf2}) | zcat - | python scripts/clean_vcf.py > {output}) 2> {log}"

    rule rename_genome_fasta:
        '''Rename the genome fasta so the other rules work with non-human genomes'''
        input: GENOME_FA
        output: KARYOTYPIC_GENOME_FA
        benchmark: "../resources/ensembl/karyotypic_order_rename.benchmark"
        log: "../resources/ensembl/karyotypic_order_rename.log"
        conda: "../envs/downloads.yaml"
        shell: "cp {input} {output} 2> {log}"

if DIVISION != "bacteria":
    rule index_ensembl_vcf:
        input: f"../resources/ensembl/{SPECIES}.ensembl.vcf"
        output: f"../resources/ensembl/{SPECIES}.ensembl.vcf.idx"
        log: f"../resources/ensembl/{SPECIES}.ensembl.vcf.idx.log"
        benchmark: f"../resources/ensembl/{SPECIES}.ensembl.vcf.idx.benchmark"
        conda: "../envs/variants.yaml"
        shell: "gatk IndexFeatureFile -I {input} 2> {log}"

rule download_chromosome_mappings:
    output: f"../resources/ChromosomeMappings/{GENOME_VERSION}_UCSC2ensembl.txt"
    params: url="https://github.com/dpryan79/ChromosomeMappings.git"
    log: "../resources/download_chromosome_mappings.log"
    benchmark: "../resources/download_chromosome_mappings.benchmark"
    conda: "../envs/downloads.yaml"
    shell: "(cd ../resources && git clone {params.url}) 2> {log}"

rule dict_fa:
    input: KARYOTYPIC_GENOME_FA
    output: f"{KARYOTYPIC_GENOME_PREFIX}.dict"
    log: f"{KARYOTYPIC_GENOME_PREFIX}.dict.log"
    benchmark: f"{KARYOTYPIC_GENOME_PREFIX}.dict.benchmark"
    conda: "../envs/variants.yaml"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output} 2> {log}"
