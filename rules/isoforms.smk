rule assemble_transcripts:
    input:
        bam="data/combined.sorted.bam",
        gff="data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".gff3"
    output: "data/combined.sorted.gtf"
    threads: 12
    log: "data/combined.sorted.gtf.log"
    shell:
        "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} 2> {log}" # strandedness: --fr for forwared or --rf for reverse

# rule filter_gtf_entries_without_strand
#     input: "data/ERR315327.sorted.gtf"
#     output: "data/ERR315327.sorted.filtered.gtf"
#     script:

rule build_gtf_sharp:
    output: "GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll"
    log: "data/GtfSharp.build.log"
    shell:
        "(cd GtfSharp && "
        "dotnet restore && "
        "dotnet build -c Release GtfSharp.sln) &> {log}"

rule filter_transcripts_add_cds:
    input:
        gtfsharp="GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll",
        gtf="data/combined.sorted.gtf",
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        refg="data/ensembl/Homo_sapiens.GRCh38.81.gff3"
    output:
        temp("data/combined.sorted.filtered.gtf"),
        "data/combined.sorted.filtered.withcds.gtf",
    shell:
        "dotnet {input.gtfsharp} -f {input.fa} -g {input.gtf} -r {input.refg}"

rule generate_snpeff_database:
    input:
        jar="SnpEff/snpEff.jar",
        gtf="data/combined.sorted.filtered.withcds.gtf",
        pfa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".pep.all.fa",
        gfa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa"
    output:
        gtf="SnpEff/data/combined.sorted.filtered.withcds.gtf/genes.gtf",
        pfa="SnpEff/data/combined.sorted.filtered.withcds.gtf/protein.fa",
        gfa="SnpEff/data/genomes/combined.sorted.filtered.withcds.gtf.fa",
    params:
        ref="combined.sorted.filtered.withcds.gtf"
    resources:
        mem_mb=16000
    log:
        "data/combined.sorted.filtered.withcds.snpeffdatabase.log"
    shell:
        "cp {input.gtf} {output.gtf} && "
        "cp {input.pfa} {output.pfa} && "
        "cp {input.gfa} {output.gfa} && "
        "echo \"\n# {params.ref}\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.genome : Human genome " + GENOME_VERSION + " using RefSeq transcripts\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.M.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.MT.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "(java -Xmx{resources.mem_mb}M -jar {input.jar} build -gtf22 -v {params.ref}) 2> {log}"
