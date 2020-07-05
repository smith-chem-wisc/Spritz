TRANSFER_MOD_DLL="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll"
REF=config["species"] + "." + config["genome"]

rule download_protein_xml:
    output:
        xml=UNIPROTXML,
        fasta=UNIPROTFASTA
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

rule generate_reference_snpeff_database:
    input:
        jar="SnpEff/snpEff.jar",
        gff3=GFF3,
        pfa="data/ensembl/{REF}.pep.all.fa",
        gfa="data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa"
    output:
        pfa="SnpEff/data/{REF}/protein.fa",
        gff3="SnpEff/data/{REF}/genes.gff",
        gfa="SnpEff/data/genomes/{REF}.fa",
        done="SnpEff/data/{REF}/done{REF}.txt"
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
        temp=directory("temporary"),
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        done="{dir}/done" + REF + "." + ENSEMBL_VERSION + ".txt",
        protxml=temp("{dir}/" + REF + "." + ENSEMBL_VERSION + ".protein.xml"),
        protxmlgz="{dir}/" + REF + "." + ENSEMBL_VERSION + ".protein.xml.gz",
        protxmlwithmods=temp("{dir}/" + REF + "." + ENSEMBL_VERSION + ".protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/" + REF + "." + ENSEMBL_VERSION + ".protein.withmods.xml.gz",
    resources: mem_mb=16000
    benchmark: "{dir}/" + REF + "." + ENSEMBL_VERSION + ".spritz.benchmark"
    log: "{dir}/" + REF + "." + ENSEMBL_VERSION + ".spritz.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {REF} && " # no isoforms, no variants
        "mv {output.protxml} {input.temp}/" + REF + "." + ENSEMBL_VERSION + ".protein.xml && "
        "dotnet {input.transfermods} -x {input.unixml} -y {input.temp}/" + REF + "." + ENSEMBL_VERSION + ".protein.xml && "
        "mv {input.temp}/" + REF + "." + ENSEMBL_VERSION + ".protein.xml {wildcards.dir} && "
        "mv {input.temp}/" + REF + "." + ENSEMBL_VERSION + ".protein.withmods.xml {wildcards.dir} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log} && touch {output.done}"

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
    resources: mem_mb=16000
    benchmark: "{dir}/combined.spritz.isoform.benchmark"
    log: "{dir}/combined.spritz.isoform.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} < /dev/null && " # isoforms, no variants
        "mv {output.protxml} {input.temp}/{params.infile} && "
        "dotnet {input.transfermods} -x {input.unixml} -y {input.temp}/{params.infile} && "
        "mv {input.temp}/{params.infile} {wildcards.dir} && "
        "mv {input.temp}/{params.outfile} {wildcards.dir} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log}"
