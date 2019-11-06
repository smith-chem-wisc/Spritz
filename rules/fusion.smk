FUSION_REF_VERSION = config["genome"] + "_gencode_v29_CTAT_lib_Mar272019.plug-n-play"

rule download_premade_fusion_indices:
    '''Get the premade STAR-Fusion indices'''
    output: "data/" + FUSION_REF_VERSION + "/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/SA"
    shell:
        "wget -O - https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/" + FUSION_REF_VERSION + ".tar.gz | "
        "tar -C data -xz"

rule unzip_for_star_fusion:
    '''Gunzip files before STAR-Fusion because it doesn't play well with gunzip commands'''
    input:
        fq1="{dir}/trimmed/{sra}.trim_1.fastq.gz" if check_sra() else "{dir}/{fq}_1.fastq.gz",
        fq2="{dir}/trimmed/{sra}.trim_2.fastq.gz" if check_sra() else "{dir}/{fq}_2.fastq.gz",
    output:
        fq1=temp("{dir}/trimmed/{sra}.trim_1.fastq") if check_sra() else "{dir}/{fq}_1.fastq",
        fq2=temp("{dir}/trimmed/{sra}.trim_2.fastq") if check_sra() else "{dir}/{fq}_2.fastq",
    threads: 2
    shell:
        "gunzip -k {input.fq1} & gunzip -k {input.fq2}"

rule rsem_star_fusion:
    '''Align with chimeric alignments and analyze coding effects with STAR-Fusion'''
    input:
        tmpdir=temp(directory("tmp")),
        genomelibsa="data/" + FUSION_REF_VERSION + "/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/SA",
        genomelibdir=directory("data/" + FUSION_REF_VERSION + "/ctat_genome_lib_build_dir"),
        fq1="{dir}/trimmed/{sra}.trim_1.fastq" if check_sra() else expand("data/{fq1}_1.fastq", fq1=config["fq1"]),
        fq2="{dir}/trimmed/{sra}.trim_2.fastq" if check_sra() else expand("data/{fq2}_2.fastq", fq2=config["fq2"]),
    output:
        "{dir}/{sra}FusionAnalysis/star-fusion.fusion_predictions.abridged.coding_effect.tsv",
        "{dir}/{sra}FusionAnalysis/Aligned.out.bam",
        "{dir}/{sra}FusionAnalysis/Chimeric.out.junction",
        temp(directory("{dir}/{sra}FusionAnalysis/star-fusion.preliminary")),
        temp(directory("{dir}/{sra}FusionAnalysis/_STARgenome")),
        temp(directory("{dir}/{sra}FusionAnalysis/_STARpass1")),
        temp(directory("{dir}/{sra}FusionAnalysis/_starF_checkpoints")),
    resources: mem_mb=50000
    threads: 12
    log: "{dir}/{sra}STARFusion.log"
    params: validation="--FusionInspector validate --denovo_reconstruct" # does realignment and trinity reconstructions to validate fusions
    shell:
        "(STAR-Fusion --examine_coding_effect {params.validation} --CPU {threads} --tmpdir {input.tmpdir} "
        " --genome_lib_dir {input.genomelibdir} --output_dir {dir}/{wildcards.sra}FusionAnalysis "
        " --left_fq {input.fq1} --right_fq {input.fq2}) &> {log}"

rule add_srr_to_output:
    input: "{dir}/{sra}FusionAnalysis/star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    output: "{dir}/{sra}star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    shell: "mv {input} {output}"

rule generate_fusion_proteins:
    '''Use coding effects to generate fusion proteins'''
    input:
        expand("{dir}/{sra}star-fusion.fusion_predictions.abridged.coding_effect.tsv", sra=config["sra"]),
        unixml=UNIPROTXML,
        transfermods=TRANSFER_MOD_DLL,
    output:
        "{dir}/FusionProteins.xml",
        "{dir}/FusionProteins.withmods.xml"
    shell:
        "dotnet {input.transfermods} -x {input.unixml} -f " +
        ",".join(expand("{dir}/{sra}star-fusion.fusion_predictions.abridged.coding_effect.tsv", sra=config["sra"]))
