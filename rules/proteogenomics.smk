UNIPROTXML="data/uniprot/" + config["species"] + ".protein.xml.gz" #"data/Homo_sapiens_202022.xml.gz"
TRANSFER_MOD_DLL="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll"
REF=config["species"] + "." + config["genome"]

rule download_protein_xml:
    output: UNIPROTXML
    shell: "python scripts/get_proteome.py && python scripts/download_xml.py | gzip -c > {output}" #fixme

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
        protxml="{dir}/combined.spritz.snpeff.protein.xml"
    output:
        protxml=temp("{dir}/combined.spritz.snpeff.protein.withmods.xml"),
        protxmlgz="{dir}/combined.spritz.snpeff.protein.withmods.xml.gz"
    log: "{dir}/combined.spritz.snpeff.protein.withmods.log"
    shell:
        "(mv {input.protxml} temporary/combined.spritz.snpeff.protein.xml && "
        "dotnet {input.transfermods} -x {input.unixml} -y temporary/combined.spritz.snpeff.protein.xml && "
        "mv temporary/* {wildcards.dir} && "
        "gzip -k {output.protxml}) &> {log}"

rule transfer_modifications_isoformvariant:
    input:
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
        protxml="{dir}/combined.spritz.isoformvariants.protein.xml"
    output:
        protxml=temp("{dir}/combined.spritz.isoformvariants.protein.withmods.xml"),
        protxmlgz="{dir}/combined.spritz.isoformvariants.protein.withmods.xml.gz"
    log: "{dir}/combined.spritz.isoformvariants.protein.withmods.log"
    shell:
        "(mv {input.protxml} temporary/combined.spritz.isoformvariants.protein.xml && "
        "dotnet {input.transfermods} -x {input.unixml} -y temporary/combined.spritz.isoformvariants.protein.xml && "
        "mv temporary/* {wildcards.dir} && "
        "gzip -k {output.protxml}) &> {log}"

rule reference_protein_xml:
    """
    Create protein XML with sequences from the reference gene model.
    """
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        dummy="{dir}/dummy.txt",
        protxml=temp("{dir}/" + config["genome"] + "." + config["snpeff"] + ".protein.xml"),
        protxmlgz="{dir}/" + config["genome"] + "." + config["snpeff"] + ".protein.xml.gz",
        protxmlwithmods=temp("{dir}/" + config["genome"] + "." + config["snpeff"] + ".protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/" + config["genome"] + "." + config["snpeff"] + ".protein.withmods.xml.gz",
    params:
        ref=config["genome"] + "." + config["snpeff"], # no isoform reconstruction
    resources:
        mem_mb=16000
    log:
        "{dir}/" + config["genome"] + "." + config["snpeff"] + ".spritz.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} && " # no isoforms, no variants
        "mv {output.protxml} temporary/GRCm38.86.protein.xml && "
        "dotnet {input.transfermods} -x {input.unixml} -y temporary/GRCm38.86.protein.xml && "
        "mv temporary/* {wildcards.dir} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log} && touch {output.dummy}"

rule custom_protein_xml:
    """
    Create protein XML with sequences from the isoform discovery gene model.
    """
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        isoform_reconstruction="SnpEff/data/combined.sorted.filtered.withcds.gtf/genes.gtf",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        protxml=temp("{dir}/combined.spritz.isoform.protein.xml"),
        protxmlgz="{dir}/combined.spritz.isoform.protein.xml.gz",
        protxmlwithmods=temp("{dir}/combined.spritz.isoform.protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/combined.spritz.isoform.protein.withmods.xml.gz"
    params:
        ref="combined.sorted.filtered.withcds.gtf" # with isoforms
    resources:
        mem_mb=16000
    log:
        "{dir}/combined.spritz.isoform.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} && " # isoforms, no variants
        "mv {output.protxml} temporary/combined.spritz.isoform.protein.xml && "
        "dotnet {input.transfermods} -x {input.unixml} -y temporary/combined.spritz.isoform.protein.xml && "
        "mv temporary/* {wildcards.dir} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log}"
