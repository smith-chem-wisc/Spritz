TRANSFER_MOD_DLL="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll"
REF=config["species"] + "." + config["genome"]

rule download_protein_xml:
    output:
        xml=UNIPROTXML,
        fasta=UNIPROTFASTA
    shell:
        "python scripts/get_proteome.py && "
        "python scripts/download_uniprot.py xml | gzip -c > {output.xml} && " #fixme
        "python scripts/download_uniprot.py fasta > {output.fasta}"

rule build_transfer_mods:
    output: TRANSFER_MOD_DLL
    log: "data/TransferUniProtModifications.build.log"
    shell:
        "(cd TransferUniProtModifications && "
        "dotnet restore && "
        "dotnet build -c Release TransferUniProtModifications.sln) &> {log}"

rule transfer_modifications_variant:
    input:
        temp=directory("temporary"),
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
        protxml="{dir}/combined.spritz.snpeff.protein.xml"
    output:
        protxml=temp("{dir}/combined.spritz.snpeff.protein.withmods.xml"),
        protxmlgz="{dir}/combined.spritz.snpeff.protein.withmods.xml.gz"
    params:
        infile="combined.spritz.snpeff.protein.xml",
        outfile="combined.spritz.snpeff.protein.withmods.xml"
    log: "{dir}/combined.spritz.snpeff.protein.withmods.log"
    shell:
        "(mv {input.protxml} {input.temp}/{params.infile} && "
        "dotnet {input.transfermods} -x {input.unixml} -y {input.temp}/{params.infile} && "
        "mv {input.temp}/{params.infile} {wildcards.dir} && "
        "mv {input.temp}/{params.outfile} {wildcards.dir} && "
        "gzip -k {output.protxml}) &> {log}"

rule transfer_modifications_isoformvariant:
    input:
        temp=directory("temporary"),
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
        protxml="{dir}/combined.spritz.isoformvariants.protein.xml"
    output:
        protxml=temp("{dir}/combined.spritz.isoformvariants.protein.withmods.xml"),
        protxmlgz="{dir}/combined.spritz.isoformvariants.protein.withmods.xml.gz"
    params:
        infile="combined.spritz.isoformvariants.protein.xml",
        outfile="combined.spritz.isoformvariants.protein.withmods.xml"
    log: "{dir}/combined.spritz.isoformvariants.protein.withmods.log"
    shell:
        "(mv {input.protxml} {input.temp}/{params.infile} && "
        "dotnet {input.transfermods} -x {input.unixml} -y {input.temp}/{params.infile} && "
        "mv {input.temp}/{params.infile} {wildcards.dir} && "
        "mv {input.temp}/{params.outfile} {wildcards.dir} && "
        "gzip -k {output.protxml}) &> {log}"

rule reference_protein_xml:
    """
    Create protein XML with sequences from the reference gene model.
    """
    input:
        "data/SnpEffDatabases.txt",
        temp=directory("temporary"),
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
        "mv {output.protxml} {input.temp}/{params.ref}.protein.xml && "
        "dotnet {input.transfermods} -x {input.unixml} -y {input.temp}/{params.ref}.protein.xml && "
        "mv {input.temp}/{params.ref}.protein.xml {wildcards.dir} && "
        "mv {input.temp}/{params.ref}.protein.withmods.xml {wildcards.dir} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log} && touch {output.dummy}"

rule custom_protein_xml:
    """
    Create protein XML with sequences from the isoform discovery gene model.
    """
    input:
        "data/SnpEffDatabases.txt",
        temp=directory("temporary"),
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        isoform_reconstruction=[
            "SnpEff/data/combined.transcripts.genome.gff3/genes.gff",
            "SnpEff/data/combined.transcripts.genome.gff3/protein.fa",
            "SnpEff/data/genomes/combined.transcripts.genome.gff3.fa",
            "SnpEff/data/combined.transcripts.genome.gff3/done.txt"],
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        protxml=temp("{dir}/combined.spritz.isoform.protein.xml"),
        protxmlgz="{dir}/combined.spritz.isoform.protein.xml.gz",
        protxmlwithmods=temp("{dir}/combined.spritz.isoform.protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/combined.spritz.isoform.protein.withmods.xml.gz"
    params:
        ref="combined.transcripts.genome.gff3", # with isoforms
        infile="combined.spritz.isoform.protein.xml",
        outfile="combined.spritz.isoform.protein.withmods.xml"
    resources:
        mem_mb=16000
    log:
        "{dir}/combined.spritz.isoform.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} < /dev/null && " # isoforms, no variants
        "mv {output.protxml} {input.temp}/{params.infile} && "
        "dotnet {input.transfermods} -x {input.unixml} -y {input.temp}/{params.infile} && "
        "mv {input.temp}/{params.infile} {wildcards.dir} && "
        "mv {input.temp}/{params.outfile} {wildcards.dir} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log}"
