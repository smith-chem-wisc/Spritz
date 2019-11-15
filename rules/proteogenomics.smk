UNIPROTXML="data/uniprot/human.protein.xml.gz" #"data/Homo_sapiens_202022.xml.gz"
TRANSFER_MOD_DLL="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll"

rule download_protein_xml:
    output: UNIPROTXML
    shell: "python scripts/download_protein_xml.py | gzip -c > {output}"

rule build_transfer_mods:
    output: TRANSFER_MOD_DLL
    log: "data/TransferUniProtModifications.build.log"
    shell:
        "(cd TransferUniProtModifications && "
        "dotnet restore && "
        "dotnet build -c Release TransferUniProtModifications.sln) &> {log}"

rule transfer_modifications_variant:
    input:
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
        protxml="data/combined.spritz.snpeff.protein.xml"
    output:
        protxml=temp("data/combined.spritz.snpeff.protein.withmods.xml"),
        protxmlgz="data/combined.spritz.snpeff.protein.withmods.xml.gz"
    log: "data/combined.spritz.snpeff.protein.withmods.log"
    shell:
        "(dotnet {input.transfermods} -x {input.unixml} -y {input.protxml} && gzip -k {output.protxml}) &> {log}" # typo

rule transfer_modifications_isoformvariant:
    input:
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
        protxml="data/combined.spritz.isoformvariants.protein.xml"
    output:
        protxml=temp("data/combined.spritz.isoformvariants.protein.withmods.xml"),
        protxmlgz="data/combined.spritz.isoformvariants.protein.withmods.xml.gz"
    log: "data/combined.spritz.isoformvariants.protein.withmods.log"
    shell:
        "(dotnet {input.transfermods} -x {input.unixml} -y {input.protxml} && gzip -k {output.protxml}) &> {log}"

rule reference_protein_xml:
    """
    Create protein XML with sequences from the reference gene model.
    """
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        protxml=temp("data/GRCh38.86.protein.xml"),
        protxmlgz="data/GRCh38.86.protein.xml.gz",
        protxmlwithmods=temp("data/GRCh38.86.protein.withmods.xml"),
        protxmlwithmodsgz="data/GRCh38.86.protein.withmods.xml.gz",
    params:
        ref=GENEMODEL_VERSION_SNPEFF, # no isoform reconstruction
    resources:
        mem_mb=16000
    log:
        "data/GRCh38.86.spritz.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} && " # no isoforms, no variants
        "dotnet {input.transfermods} -x {input.unixml} -y {output.protxml} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log}"

rule custom_protein_xml:
    """
    Create protein XML with sequences from the isoform discovery gene model.
    """
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        isoform_reconstruction="SnpEff/data/combined.sorted.filtered.withcds.gtf/genes.gtf",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        protxml=temp("data/combined.spritz.isoform.protein.xml"),
        protxmlgz="data/combined.spritz.isoform.protein.xml.gz",
        protxmlwithmods=temp("data/combined.spritz.isoform.protein.withmods.xml"),
        protxmlwithmodsgz="data/combined.spritz.isoform.protein.withmods.xml.gz"
    params:
        ref="combined.sorted.filtered.withcds.gtf" # with isoforms
    resources:
        mem_mb=16000
    log:
        "data/combined.spritz.isoform.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} && " # isoforms, no variants
        "dotnet {input.transfermods} -x {input.unixml} -y {output.protxml} &&"
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log}"
