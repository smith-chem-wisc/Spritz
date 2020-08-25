GATK_MEM=16000 # MB
GATK_JAVA=f"--java-options \"-Xmx{GATK_MEM}M -Dsamjdk.compression_level=9\""

rule download_snpeff:
    output: "SnpEff/snpEff.config", "SnpEff/snpEff.jar"
    log: "data/SnpEffInstall.log"
    shell:
        """
        (git clone --depth=1 https://github.com/smith-chem-wisc/SnpEff
        cd SnpEff
        mvn install:install-file -Dfile=lib/antlr-4.5.1-complete.jar -DgroupId=org.antlr -DartifactId=antlr -Dversion=4.5.1 -Dpackaging=jar
        mvn install:install-file -Dfile=lib/biojava3-core-3.0.7.jar -DgroupId=org.biojava -DartifactId=biojava3-core -Dversion=3.0.7 -Dpackaging=jar
        mvn install:install-file -Dfile=lib/biojava3-structure-3.0.7.jar -DgroupId=org.biojava -DartifactId=biojava3-structure -Dversion=3.0.7 -Dpackaging=jar
        export VERSION=4.3
        export VERSION_UND=`echo $VERSION | tr '.' '_'`
        mvn clean compile assembly:assembly
        mvn install:install-file -Dfile=target/SnpEff-$VERSION.jar -DgroupId=org.snpeff -DartifactId=SnpEff -Dversion=$VERSION -Dpackaging=jar -DgeneratePom=true --quiet
        cp target/SnpEff-$VERSION-jar-with-dependencies.jar snpEff.jar
        cd ..) &> {log}
        """

rule index_fa:
    input: "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa.fai"
    shell: "samtools faidx data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa"

rule dict_fa:
    input: "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.dict"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"

rule tmpdir:
    output: temp(directory("tmp{q}"))
    shell: "mkdir tmp{wildcards.q}"

rule hisat2_groupmark_bam:
    input:
        sorted="data/{q}/combined.sorted.bam",
        tmp=directory("tmp{q}")
    output:
        grouped=temp("data/{q}/combined.sorted.grouped.bam"),
        groupedidx=temp("data/{q}/combined.sorted.grouped.bam.bai"),
        marked="data/{q}/combined.sorted.grouped.marked.bam",
        markedidx="data/{q}/combined.sorted.grouped.marked.bam.bai",
        metrics="data/{q}/combined.sorted.grouped.marked.metrics"
    resources: mem_mb=GATK_MEM
    log: "data/{q}/combined.sorted.grouped.marked.log"
    shell:
        "(gatk {GATK_JAVA} AddOrReplaceReadGroups -PU platform  -PL illumina -SM sample -LB library -I {input.sorted} -O {output.grouped} -SO coordinate --TMP_DIR {input.tmp} && "
        "samtools index {output.grouped} && "
        "gatk {GATK_JAVA} MarkDuplicates -I {output.grouped} -O {output.marked} -M {output.metrics} --TMP_DIR {input.tmp} -AS true &&"
        "samtools index {output.marked}) &> {log}"

# Checks if quality encoding is correct, and then splits n cigar reads
rule split_n_cigar_reads:
    input:
        bam="data/{q}/combined.sorted.grouped.marked.bam",
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        fai="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa.fai",
        fadict="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.dict",
        tmp=directory("tmp{q}")
    output:
        fixed=temp("data/{q}/combined.fixedQuals.bam"),
        split=temp("data/{q}/combined.sorted.grouped.marked.split.bam"),
        splitidx=temp("data/{q}/combined.sorted.grouped.marked.split.bam.bai")
    resources: mem_mb=GATK_MEM
    log: "data/{q}/combined.sorted.grouped.marked.split.log"
    shell:
        "(gatk {GATK_JAVA} FixMisencodedBaseQualityReads -I {input.bam} -O {output.fixed} && "
        "gatk {GATK_JAVA} SplitNCigarReads -R {input.fa} -I {output.fixed} -O {output.split} --tmp-dir {input.tmp} || " # fix and split
        "gatk {GATK_JAVA} SplitNCigarReads -R {input.fa} -I {input.bam} -O {output.split} --tmp-dir {input.tmp}; " # or just split
        "samtools index {output.split}) &> {log}" # always index

rule base_recalibration:
    input:
        knownsites="data/ensembl/common_all_20170710.ensembl.vcf",
        knownsitesidx="data/ensembl/common_all_20170710.ensembl.vcf.idx",
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        bam="data/{q}/combined.sorted.grouped.marked.split.bam",
        tmp=directory("tmp{q}")
    output:
        recaltable=temp("data/{q}/combined.sorted.grouped.marked.split.recaltable"),
        recalbam=temp("data/{q}/combined.sorted.grouped.marked.split.recal.bam")
    resources: mem_mb=GATK_MEM
    log: "data/{q}/combined.sorted.grouped.marked.split.recal.log"
    shell:
        "(gatk {GATK_JAVA} BaseRecalibrator -R {input.fa} -I {input.bam} --known-sites {input.knownsites} -O {output.recaltable} --tmp-dir {input.tmp} && "
        "gatk {GATK_JAVA} ApplyBQSR -R {input.fa} -I {input.bam} --bqsr-recal-file {output.recaltable} -O {output.recalbam} --tmp-dir {input.tmp} && "
        "samtools index {output.recalbam}) &> {log}"

rule call_gvcf_varaints:
    input:
        knownsites="data/ensembl/common_all_20170710.ensembl.vcf",
        knownsitesidx="data/ensembl/common_all_20170710.ensembl.vcf.idx",
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        bam="data/{q}/combined.sorted.grouped.marked.split.recal.bam",
        tmp=directory("tmp{q}")
    output: temp("data/{q}/combined.sorted.grouped.marked.split.recal.g.vcf.gz"),
    threads: 8
        # HaplotypeCaller is only fairly efficient with threading;
        # ~14000 regions/min with 24 threads,
        # and ~13000 regions/min with 8 threads,
        # so going with 8 threads max here
    resources: mem_mb=GATK_MEM
    log: "data/{q}/combined.sorted.grouped.marked.split.recal.g.log"
    shell:
        "(gatk {GATK_JAVA} HaplotypeCaller"
        " --native-pair-hmm-threads {threads}"
        " -R {input.fa} -I {input.bam}"
        " --min-base-quality-score 20 --dont-use-soft-clipped-bases true"
        " --dbsnp {input.knownsites} -O {output} --tmp-dir {input.tmp}"
        " -ERC GVCF --max-mnp-distance 3 &&"
        "gatk IndexFeatureFile -F {output}) &> {log}"

rule call_vcf_variants:
    input:
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        gvcf="data/{q}/combined.sorted.grouped.marked.split.recal.g.vcf.gz",
        tmp=directory("tmp{q}")
    output: "data/{q}/combined.sorted.grouped.marked.split.recal.g.gt.vcf" # renamed in next rule
    resources: mem_mb=GATK_MEM
    log: "data/{q}/combined.sorted.grouped.marked.split.recal.g.gt.log"
    shell:
        "(gatk {GATK_JAVA} GenotypeGVCFs -R {input.fa} -V {input.gvcf} -O {output} --tmp-dir {input.tmp} && "
        "gatk IndexFeatureFile -F {output}) &> {log}"

rule final_vcf_naming:
    input: "data/{q}/combined.sorted.grouped.marked.split.recal.g.gt.vcf"
    output: "data/{q}/combined.spritz.vcf"
    shell: "mv {input} {output}"

rule snpeff_database_setup:
    input:
        # dir="SnpEff/data",
        jar="SnpEff/snpEff.jar",
        config="SnpEff/snpEff.config"
    output: "data/SnpEffDatabases.txt"
    params: ref=GENEMODEL_VERSION_SNPEFF
    resources: mem_mb=16000
    shell:
        "java -Xmx{resources.mem_mb}M -jar {input.jar} databases > {output} && "
        "echo \"\n# {params.ref}\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.genome : Human genome " + GENOME_VERSION + " using RefSeq transcripts\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.M.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.MT.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config"

rule variant_annotation_ref:
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        vcf="data/{q}/combined.spritz.vcf",
    output:
        ann="data/{q}/combined.spritz.snpeff.vcf",
        html="data/{q}/combined.spritz.snpeff.html",
        genesummary="data/{q}/combined.spritz.snpeff.genes.txt",
        protfa="data/{q}/combined.spritz.snpeff.protein.fasta",
        protxml="data/{q}/combined.spritz.snpeff.protein.xml"
    params: ref=GENEMODEL_VERSION_SNPEFF, # no isoform reconstruction
    resources: mem_mb=GATK_MEM
    log: "data/{q}/combined.spritz.snpeff.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml} "
        " {params.ref} {input.vcf}" # no isoforms, with variants
        " > {output.ann}) 2> {log}"
