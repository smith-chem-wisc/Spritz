TRANSFER_MOD_DLL="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll"
REF=config["species"] + "." + config["genome"]

rule download_protein_xml:
    output:
        xml=UNIPROTXML,
        fasta=UNIPROTFASTA,
    log: UNIPROTXML + ".log"
    shell:
        "(python scripts/get_proteome.py && "
        "python scripts/download_uniprot.py xml | gzip -c > {output.xml} && " #fixme
        "python scripts/download_uniprot.py fasta > {output.fasta}) &> {log}"

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
        protxml="{dir}/variants/combined.spritz.snpeff.protein.xml",
    output:
        protfastawithdecoys="{dir}/variants/combined.spritz.snpeff.protein.withdecoys.fasta",
        protxmlgz="{dir}/variants/combined.spritz.snpeff.protein.xml.gz",
        protxmlwithmods=temp("{dir}/variants/combined.spritz.snpeff.protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/variants/combined.spritz.snpeff.protein.withmods.xml.gz",
    log: "{dir}/combined.spritz.snpeff.protein.withmods.log"
    shell:
        "(dotnet {input.transfermods} -x {input.unixml} -y {input.protxml} && "
        "gzip -k {input.protxml} {output.protxmlwithmods}) &> {log}"

rule transfer_modifications_isoformvariant:
    input:
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
        protxml="{dir}/variants/combined.spritz.isoformvariants.protein.xml",
    output:
        protfastawithdecoys="{dir}/variants/combined.spritz.isoformvariants.protein.withdecoys.fasta",
        protxmlgz="{dir}/variants/combined.spritz.isoformvariants.protein.xml.gz",
        protxmlwithmods=temp("{dir}/variants/combined.spritz.isoformvariants.protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/variants/combined.spritz.isoformvariants.protein.withmods.xml.gz",
    log: "{dir}/combined.spritz.isoformvariants.protein.withmods.log"
    shell:
        "(dotnet {input.transfermods} -x {input.unixml} -y {input.protxml} && "
        "gzip -k {output.protxmlwithmods} {input.protxml}) &> {log}"

rule generate_reference_snpeff_database:
    input:
        jar="SnpEff/snpEff.jar",
        gff3=GFF3,
        pfa="data/ensembl/{REF}.pep.all.fa",
        gfa="data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
    output:
        pfa="SnpEff/data/{REF}/protein.fa",
        gff3="SnpEff/data/{REF}/genes.gff",
        gfa="SnpEff/data/genomes/{REF}.fa",
        done="SnpEff/data/{REF}/done{REF}.txt",
    resources: mem_mb=16000
    benchmark: "SnpEff/data/{REF}/snpeffdatabase.benchmark"
    log: "SnpEff/data/{REF}/snpeffdatabase.log"
    shell:
        "cp {input.gff3} {output.gff3} && "
        "cp {input.pfa} {output.pfa} && "
        "cp {input.gfa} {output.gfa} && "
        "echo \"\n# {REF}\" >> SnpEff/snpEff.config && "
        "echo \"{REF}.genome : Human genome " + GENOME_VERSION + " using RefSeq transcripts\" >> SnpEff/snpEff.config && "
        "echo \"{REF}.reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> SnpEff/snpEff.config && "
        "echo \"\t{REF}.M.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "echo \"\t{REF}.MT.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "(java -Xmx{resources.mem_mb}M -jar {input.jar} build -gff3 -v {REF}) &> {log} && touch {output.done}"

rule reference_protein_xml:
    """
    Create protein XML with sequences from the reference gene model.
    """
    input:
        "SnpEff/data/" + REF + "/done" + REF + ".txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        done="{dir}/variants/done" + REF + "." + ENSEMBL_VERSION + ".txt",
        protxml=temp("{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.xml"),
        protxmlgz="{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.xml.gz",
        protfa="{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.fasta",
        protwithdecoysfa="{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.withdecoys.fasta",
        protxmlwithmods=temp("{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.withmods.xml.gz",
    resources: mem_mb=16000
    benchmark: "{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".spritz.benchmark"
    log: "{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".spritz.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {REF} && " # no isoforms, no variants
        "dotnet {input.transfermods} -x {input.unixml} -y {output.protxml} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log} && touch {output.done}"

rule custom_protein_xml:
    """
    Create protein XML with sequences from the isoform discovery gene model.
    """
    input:
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
        protxml=temp("{dir}/isoforms/combined.spritz.isoform.protein.xml"),
        protwithdecoysfa="{dir}/isoforms/combined.spritz.isoform.protein.withdecoys.fasta",
        protxmlgz="{dir}/isoforms/combined.spritz.isoform.protein.xml.gz",
        protxmlwithmods=temp("{dir}/isoforms/combined.spritz.isoform.protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/isoforms/combined.spritz.isoform.protein.withmods.xml.gz",
        protfa="{dir}/isoforms/combined.spritz.isoform.protein.fasta",
    params:
        ref="combined.transcripts.genome.gff3", # with isoforms
    resources: mem_mb=16000
    benchmark: "{dir}/isoforms/combined.spritz.isoform.benchmark"
    log: "{dir}/isoforms/combined.spritz.isoform.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} < /dev/null && " # isoforms, no variants
        "dotnet {input.transfermods} -x {input.unixml} -y {output.protxml} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log}"
