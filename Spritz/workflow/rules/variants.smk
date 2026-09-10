GATK_MEM=24000 # MB
GATK_JAVA=f"--java-options \"-Xmx{GATK_MEM}M -Dsamjdk.compression_level=9\""

rule download_snpeff:
    '''Download and unpack custom SnpEff for annotating variants'''
    output:
        "../resources/SnpEff/snpEff.config",
        "../resources/SnpEff/snpEff.jar",
        filename=temp("../resources/SnpEff_4.3_SmithChemWisc_v2.zip")
    params:
        url="https://github.com/smith-chem-wisc/SnpEff/releases/download/4.3_SCW1/SnpEff_4.3_SmithChemWisc_v2.zip"
    log: "../resources/SnpEffInstall.log"
    conda: "../envs/downloads.yaml"
    shell:
        "(cd ../resources/ && "
        "wget {params.url} && "
        "unzip {output.filename} -d SnpEff) &> {log}"

rule index_fa:
    '''Index genome FASTA file'''
    input: KARYOTYPIC_GENOME_FA
    output: f"{KARYOTYPIC_GENOME_PREFIX}.fa.fai"
    log: f"{KARYOTYPIC_GENOME_PREFIX}.fa.faindex.log"
    conda: "../envs/variants.yaml"
    shell: "samtools faidx {input}"

rule variant_tmpdir:
    output: temp(directory("../resources/tmp")),
    log: "../resources/tmpdir.log"
    conda: "../envs/variants.yaml"
    shell: "mkdir {output} 2> {log}"

rule hisat2_group:
    '''Add read groups to sorted BAM file'''
    input:
        sorted="{dir}/align/combined.sorted.bam",
        tmp="../resources/tmp"
    output:
        grouped=temp("{dir}/variants/combined.sorted.grouped.bam"),
        groupedidx=temp("{dir}/variants/combined.sorted.grouped.bam.bai")
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} AddOrReplaceReadGroups"
        " -PU platform -PL illumina -SM sample -LB library"
        " -I {input.sorted} -O {output.grouped} --TMP_DIR {input.tmp} && "
        "samtools index {output.grouped}) &> {log}"

rule hisat2_mark:
    '''Mark duplicates in sorted BAM file'''
    input:
        grouped="{dir}/variants/combined.sorted.grouped.bam",
        tmp="../resources/tmp"
    output:
        marked="{dir}/variants/combined.sorted.grouped.marked.bam",
        markedidx="{dir}/variants/combined.sorted.grouped.marked.bam.bai",
        metrics="{dir}/variants/combined.sorted.grouped.marked.metrics"
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} MarkDuplicates"
        " -I {input.grouped} -O {output.marked} -M {output.metrics}"
        " --ASSUME_SORT_ORDER coordinate --TMP_DIR {input.tmp} && "
        "samtools index {output.marked}) &> {log}"

# Checks if quality encoding is correct, and then splits n cigar reads
rule split_n_cigar_reads:
    '''Check quality scores and split Ns in the cigar reads'''
    input:
        bam="{dir}/variants/combined.sorted.grouped.marked.bam",
        fa=KARYOTYPIC_GENOME_FA,
        fai=f"{KARYOTYPIC_GENOME_PREFIX}.fa.fai",
        fadict=f"{KARYOTYPIC_GENOME_PREFIX}.dict",
        tmp="../resources/tmp"
    output:
        fixed=temp("{dir}/variants/combined.fixedQuals.bam"),
        split=temp("{dir}/variants/combined.sorted.grouped.marked.split.bam"),
        splitidx=temp("{dir}/variants/combined.sorted.grouped.marked.split.bam.bai")
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} FixMisencodedBaseQualityReads -I {input.bam} -O {output.fixed} && "
        "gatk {params.gatk_java} SplitNCigarReads -R {input.fa} -I {output.fixed} -O {output.split} --tmp-dir {input.tmp} || " # fix and split
        "gatk {params.gatk_java} SplitNCigarReads -R {input.fa} -I {input.bam} -O {output.split} --tmp-dir {input.tmp}; " # or just split
        "samtools index {output.split}) &> {log}" # always index

if KNOWN_SITES == "bootstrap":
    # GATK's advice when no known set exists: call once on unrecalibrated data, keep the calls you
    # trust most, recalibrate against those, then call for real. Issue #186.
    #
    # One round, not iterated to convergence. Snakemake builds a static DAG, so repeating would mean
    # either unrolling a fixed number of passes or introducing checkpoints, and a single round is
    # what GATK's own RNA-seq workflow does. Convergence would be a later change.
    rule bootstrap_call_variants:
        '''First-pass calling on unrecalibrated alignments, only to obtain known sites'''
        input:
            fa=KARYOTYPIC_GENOME_FA,
            fai=f"{KARYOTYPIC_GENOME_PREFIX}.fa.fai",
            fadict=f"{KARYOTYPIC_GENOME_PREFIX}.dict",
            bam="{dir}/variants/combined.sorted.grouped.marked.split.bam",
            tmp="../resources/tmp"
        output: temp("{dir}/variants/bootstrap.g.vcf.gz")
        threads: 8
        params: gatk_java=GATK_JAVA
        resources: mem_mb=GATK_MEM
        log: "{dir}/variants/bootstrap.g.log"
        benchmark: "{dir}/variants/bootstrap.g.benchmark"
        conda: "../envs/variants.yaml"
        # No --dbsnp: that only stamps rsIDs from a known set, which is the thing missing here.
        shell:
            "(gatk {params.gatk_java} HaplotypeCaller"
            " --native-pair-hmm-threads {threads}"
            " -R {input.fa} -I {input.bam}"
            " --min-base-quality-score 20 --dont-use-soft-clipped-bases true"
            " -O {output} --tmp-dir {input.tmp}"
            " -ERC GVCF --max-mnp-distance 3 && "
            "gatk IndexFeatureFile -I {output}) &> {log}"

    rule bootstrap_genotype:
        '''Genotype the first pass so it can be filtered'''
        input:
            fa=KARYOTYPIC_GENOME_FA,
            gvcf="{dir}/variants/bootstrap.g.vcf.gz",
            tmp="../resources/tmp"
        output: temp("{dir}/variants/bootstrap.gt.vcf")
        params: gatk_java=GATK_JAVA
        resources: mem_mb=GATK_MEM
        log: "{dir}/variants/bootstrap.gt.log"
        benchmark: "{dir}/variants/bootstrap.gt.benchmark"
        conda: "../envs/variants.yaml"
        shell:
            "(gatk {params.gatk_java} GenotypeGVCFs"
            " -R {input.fa} -V {input.gvcf} -O {output} --tmp-dir {input.tmp} && "
            "gatk IndexFeatureFile -I {output}) &> {log}"

    rule bootstrap_known_sites:
        '''Keep only the first-pass SNPs worth recalibrating against'''
        input:
            fa=KARYOTYPIC_GENOME_FA,
            vcf="{dir}/variants/bootstrap.gt.vcf",
            tmp="../resources/tmp"
        output:
            filtered=temp("{dir}/variants/bootstrap.filtered.vcf"),
            vcf="{dir}/variants/bootstrap.knownsites.vcf",
            idx="{dir}/variants/bootstrap.knownsites.vcf.idx",
        params: gatk_java=GATK_JAVA
        resources: mem_mb=GATK_MEM
        log: "{dir}/variants/bootstrap.knownsites.log"
        benchmark: "{dir}/variants/bootstrap.knownsites.benchmark"
        conda: "../envs/variants.yaml"
        # GATK's RNA-seq hard filters, which are not the DNA ones: a 35-base window allowing 3
        # clustered SNPs, FS > 30 rather than 60, and QD < 2. Spritz aligns RNA-seq with hisat2 and
        # has already run SplitNCigarReads, so the RNA-seq thresholds are the applicable ones.
        #
        # Then indels and anything filtered are dropped, and QUAL >= 30 keeps this to sites worth
        # calling known. Recalibration treats every non-known mismatch as an error, so a permissive
        # set here is worse than a small one.
        shell:
            "(gatk {params.gatk_java} VariantFiltration"
            " -R {input.fa} -V {input.vcf} -O {output.filtered}"
            " --window 35 --cluster 3"
            " --filter-name FS -filter \"FS > 30.0\""
            " --filter-name QD -filter \"QD < 2.0\" --tmp-dir {input.tmp} && "
            "gatk {params.gatk_java} SelectVariants"
            " -R {input.fa} -V {output.filtered} -O {output.vcf}"
            " --select-type-to-include SNP --exclude-filtered true"
            " -select \"QUAL >= 30.0\" --tmp-dir {input.tmp} && "
            "gatk IndexFeatureFile -I {output.vcf}) &> {log}"

rule base_recalibration:
    '''Generate recalibration table and recalibrate BAM file'''
    input:
        knownsites=KNOWN_SITES_VCF,
        knownsitesidx=KNOWN_SITES_VCF_IDX,
        fa=KARYOTYPIC_GENOME_FA,
        bam="{dir}/variants/combined.sorted.grouped.marked.split.bam",
        tmp="../resources/tmp"
    output:
        recaltable=temp("{dir}/variants/combined.sorted.grouped.marked.split.recaltable"),
        recalbam=temp("{dir}/variants/combined.sorted.grouped.marked.split.recal.bam")
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} BaseRecalibrator -R {input.fa} -I {input.bam}"
        " --known-sites {input.knownsites} -O {output.recaltable} --tmp-dir {input.tmp} && "
        "gatk {params.gatk_java} ApplyBQSR -R {input.fa} -I {input.bam}"
        " --bqsr-recal-file {output.recaltable} -O {output.recalbam} --tmp-dir {input.tmp} && "
        "samtools index {output.recalbam}) &> {log}"

rule call_gvcf_varaints:
    '''Create genome VCF file'''
    input:
        # Declared in both modes so the recalibrated BAM and this call see the same set, but only
        # passed to --dbsnp in ensembl mode; see params.dbsnp.
        knownsites=KNOWN_SITES_VCF,
        knownsitesidx=KNOWN_SITES_VCF_IDX,
        fa=KARYOTYPIC_GENOME_FA,
        bam="{dir}/variants/combined.sorted.grouped.marked.split.recal.bam",
        tmp="../resources/tmp"
    output: temp("{dir}/variants/combined.sorted.grouped.marked.split.recal.g.vcf.gz"),
    threads: 8
        # HaplotypeCaller is only fairly efficient with threading;
        # ~14000 regions/min with 24 threads,
        # and ~13000 regions/min with 8 threads,
        # so going with 8 threads max here
    params:
        gatk_java=GATK_JAVA,
        # rsIDs only. A bootstrap set is this run's own first-pass calls, so stamping them as known
        # identifiers would invent provenance that does not exist.
        dbsnp=(lambda w, input: f"--dbsnp {input.knownsites}") if KNOWN_SITES == "ensembl" else ""
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} HaplotypeCaller"
        " --native-pair-hmm-threads {threads}"
        " -R {input.fa} -I {input.bam}"
        " --min-base-quality-score 20 --dont-use-soft-clipped-bases true"
        " {params.dbsnp} -O {output} --tmp-dir {input.tmp}"
        " -ERC GVCF --max-mnp-distance 3 && "
        "gatk IndexFeatureFile -I {output}) &> {log}"

rule call_vcf_variants:
    '''Genotype the gVCF for the combined dataset to make VCF'''
    input:
        fa=KARYOTYPIC_GENOME_FA,
        gvcf="{dir}/variants/combined.sorted.grouped.marked.split.recal.g.vcf.gz",
        tmp="../resources/tmp"
    output: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.vcf" # renamed in next rule
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} GenotypeGVCFs"
        " -R {input.fa} -V {input.gvcf} -O {output} --tmp-dir {input.tmp} && "
        "gatk IndexFeatureFile -I {output}) &> {log}"

# Both rules below produce combined.spritz.vcf, the point every downstream rule reads the variant set
# from, so they are mutually exclusive rather than alternatives snakemake picks between: declaring one
# output twice is an AmbiguousRuleException, and resolving it with ruleorder would leave the losing
# branch's whole input chain still in the DAG.
if check('vcf'):
    # Only defined when there is something to merge. With a single VCF the file the user supplied is
    # the input to stage_user_vcf directly, so nothing rewrites it and the bcftools env is never even
    # built - the one-file case behaves exactly as it did before merging existed.
    if len(config['vcf']) > 1:
        rule merge_user_vcfs:
            '''Combine per-sample VCFs into the one multi-sample VCF the workflow reads'''
            input:
                # Relative to the analysis directory, like fq and fq_se: only that directory and
                # resources/ are bind-mounted into the container, so a host path would not resolve.
                vcfs=lambda w: [posixpath.join(w.dir, v) for v in config['vcf']],
            output: temp("{dir}/variants/user.supplied.vcf")
            log: "{dir}/variants/merge_user_vcfs.log"
            benchmark: "{dir}/variants/merge_user_vcfs.benchmark"
            conda: "../envs/variants.yaml"
            # bcftools rather than GATK: MergeVcfs requires every input to carry the same sample set,
            # so it cannot combine one-sample-per-file VCFs, which is the case this exists for.
            #
            # Sorted on the way in because `index -t` fails on an unsorted VCF, and a caller that
            # emitted one is a bad error message rather than a bad input. Merged in the order given
            # rather than by globbing the temporary directory, which would order 10 before 2.
            #
            # --force-samples renames a collision instead of aborting. Callers routinely emit a
            # placeholder sample name, so two files both naming their sample SAMPLE is ordinary.
            shell:
                "(tmp=$(mktemp -d) && i=0 && sorted= && "
                "for v in {input.vcfs}; do i=$((i+1)); "
                "bcftools sort -Oz -o \"$tmp/$i.vcf.gz\" \"$v\" && "
                "bcftools index -t \"$tmp/$i.vcf.gz\" && "
                "sorted=\"$sorted $tmp/$i.vcf.gz\"; done && "
                "bcftools merge --force-samples -Ov -o {output} $sorted && "
                "rm -rf \"$tmp\") &> {log}"

        USER_SUPPLIED_VCF = "{dir}/variants/user.supplied.vcf"
    else:
        USER_SUPPLIED_VCF = lambda w: posixpath.join(w.dir, config['vcf'][0])

    rule stage_user_vcf:
        '''Use a VCF the user called elsewhere, skipping alignment and GATK'''
        input:
            vcf=USER_SUPPLIED_VCF,
            fai=f"{KARYOTYPIC_GENOME_PREFIX}.fa.fai",
        output: "{dir}/variants/combined.spritz.vcf"
        log: "{dir}/variants/stage_user_vcf.log"
        benchmark: "{dir}/variants/stage_user_vcf.benchmark"
        conda: "../envs/default.yaml"
        shell: "python scripts/stage_user_vcf.py {input.vcf} {input.fai} > {output} 2> {log}"

else:
    rule final_vcf_naming:
        '''Rename VCF to shorter filename'''
        input: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.vcf"
        output: "{dir}/variants/combined.spritz.vcf"
        log: "{dir}/variants/final_vcf_naming.log"
        conda: "../envs/variants.yaml"
        shell: "mv {input} {output} 2> {log}"

rule variant_annotation_ref:
    '''Generate proteome FASTA and XML for reference database'''
    input:
        f"../resources/SnpEff/data/{REF}/done{REF}.txt",
        snpeff="../resources/SnpEff/snpEff.jar",
        fa=KARYOTYPIC_GENOME_FA,
        vcf="{dir}/variants/combined.spritz.vcf",
    output:
        ann="{dir}/variants/combined.spritz.snpeff.vcf",
        html="{dir}/variants/combined.spritz.snpeff.html",
        genesummary="{dir}/variants/combined.spritz.snpeff.genes.txt",
        protfa="{dir}/variants/combined.spritz.snpeff.protein.fasta",
        protxml="{dir}/variants/combined.spritz.snpeff.protein.xml"
    params: ref=REF, # no isoform reconstruction
    resources: mem_mb=16000
    log: "{dir}/variants/combined.spritz.snpeff.log"
    benchmark: "{dir}/variants/combined.spritz.snpeff.benchmark"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml} "
        " {params.ref} {input.vcf}" # no isoforms, with variants
        " > {output.ann}) 2> {log}"

rule variant_annotation_custom:
    input:
        snpeff="../resources/SnpEff/snpEff.jar",
        fa=KARYOTYPIC_GENOME_FA,
        vcf="{dir}/variants/combined.spritz.vcf",
        isoform_reconstruction=[
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/genes.gff",
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/protein.fa",
            "../resources/SnpEff/data/genomes/combined.transcripts.genome.gff3.fa",
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/done.txt"],
    output:
        ann="{dir}/variants/combined.spritz.isoformvariants.vcf",
        html="{dir}/variants/combined.spritz.isoformvariants.html",
        genesummary="{dir}/variants/combined.spritz.isoformvariants.genes.txt",
        protfa="{dir}/variants/combined.spritz.isoformvariants.protein.fasta",
        protxml=temp("{dir}/variants/combined.spritz.isoformvariants.protein.xml"),
    params: ref="combined.transcripts.genome.gff3" # with isoforms
    resources: mem_mb=GATK_MEM
    log: "{dir}/variants/combined.spritz.isoformvariants.log"
    benchmark: "{dir}/variants/combined.spritz.isoformvariants.benchmark"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml}"
        " {params.ref} {input.vcf}" # with isoforms and variants
        " > {output.ann}) 2> {log}"

rule finish_variants:
    '''Copy final output files from variant workflow to main directory'''
    input:
        ann="{dir}/variants/combined.spritz.snpeff.vcf",
        protfa="{dir}/variants/combined.spritz.snpeff.protein.fasta",
        protwithdecoysfa="{dir}/variants/combined.spritz.snpeff.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/variants/combined.spritz.snpeff.protein.withmods.xml.gz",
        refprotfa=posixpath.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.fasta"),
        refprotwithdecoysfa=posixpath.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.withdecoys.fasta"),
        refprotwithmodsxml=posixpath.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.withmods.xml.gz"),
        protfragpipefa="{dir}/variants/combined.spritz.snpeff.protein.fragpipe.fasta",
        protwithdecoysfragpipefa="{dir}/variants/combined.spritz.snpeff.protein.withdecoys.fragpipe.fasta",
        accname="{dir}/variants/combined.spritz.snpeff.protein.accname.tsv",
        vardesc="{dir}/variants/combined.spritz.snpeff.protein.vardesc.tsv",
        refprotfragpipefa=posixpath.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.fragpipe.fasta"),
        refprotwithdecoysfragpipefa=posixpath.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.withdecoys.fragpipe.fasta"),
        refaccname=posixpath.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.accname.tsv"),
        refvardesc=posixpath.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.vardesc.tsv"),
    output:
        ann="{dir}/final/combined.spritz.snpeff.vcf",
        protfa="{dir}/final/combined.spritz.snpeff.protein.fasta",
        protwithdecoysfa="{dir}/final/combined.spritz.snpeff.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/final/combined.spritz.snpeff.protein.withmods.xml.gz",
        refprotfa=posixpath.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.fasta"),
        refprotwithdecoysfa=posixpath.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.withdecoys.fasta"),
        refprotwithmodsxml=posixpath.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.withmods.xml.gz"),
        protfragpipefa="{dir}/final/combined.spritz.snpeff.protein.fragpipe.fasta",
        protwithdecoysfragpipefa="{dir}/final/combined.spritz.snpeff.protein.withdecoys.fragpipe.fasta",
        accname="{dir}/final/combined.spritz.snpeff.protein.accname.tsv",
        vardesc="{dir}/final/combined.spritz.snpeff.protein.vardesc.tsv",
        refprotfragpipefa=posixpath.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.fragpipe.fasta"),
        refprotwithdecoysfragpipefa=posixpath.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.withdecoys.fragpipe.fasta"),
        refaccname=posixpath.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.accname.tsv"),
        refvardesc=posixpath.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.vardesc.tsv"),
    log: "{dir}/variants/finish_isoform_variants.log"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "cp {input.ann} {input.protfa} {input.protwithdecoysfa} {input.protxmlwithmodsgz}"
        " {input.refprotfa} {input.refprotwithdecoysfa} {input.refprotwithmodsxml}"
        " {input.protfragpipefa} {input.protwithdecoysfragpipefa} {input.accname} {input.vardesc}"
        " {input.refprotfragpipefa} {input.refprotwithdecoysfragpipefa} {input.refaccname} {input.refvardesc}"
        " {wildcards.dir}/final 2> {log}"

rule finish_isoform_variants:
    '''Copy final output files from isoform-variant workflow to main directory'''
    input:
        ann="{dir}/variants/combined.spritz.isoformvariants.vcf",
        protfa="{dir}/variants/combined.spritz.isoformvariants.protein.fasta",
        protwithdecoysfa="{dir}/variants/combined.spritz.isoformvariants.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/variants/combined.spritz.isoformvariants.protein.withmods.xml.gz",
        protfragpipefa="{dir}/variants/combined.spritz.isoformvariants.protein.fragpipe.fasta",
        protwithdecoysfragpipefa="{dir}/variants/combined.spritz.isoformvariants.protein.withdecoys.fragpipe.fasta",
        accname="{dir}/variants/combined.spritz.isoformvariants.protein.accname.tsv",
        vardesc="{dir}/variants/combined.spritz.isoformvariants.protein.vardesc.tsv",
    output:
        ann="{dir}/final/combined.spritz.isoformvariants.vcf",
        protfa="{dir}/final/combined.spritz.isoformvariants.protein.fasta",
        protwithdecoysfa="{dir}/final/combined.spritz.isoformvariants.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/final/combined.spritz.isoformvariants.protein.withmods.xml.gz",
        protfragpipefa="{dir}/final/combined.spritz.isoformvariants.protein.fragpipe.fasta",
        protwithdecoysfragpipefa="{dir}/final/combined.spritz.isoformvariants.protein.withdecoys.fragpipe.fasta",
        accname="{dir}/final/combined.spritz.isoformvariants.protein.accname.tsv",
        vardesc="{dir}/final/combined.spritz.isoformvariants.protein.vardesc.tsv",
    log: "{dir}/variants/finish_isoform_variants.log"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "cp {input.ann} {input.protfa} {input.protwithdecoysfa} {input.protxmlwithmodsgz}"
        " {input.protfragpipefa} {input.protwithdecoysfragpipefa} {input.accname} {input.vardesc} {wildcards.dir}/final 2> {log}"
