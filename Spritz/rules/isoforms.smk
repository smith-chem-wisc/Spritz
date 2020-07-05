import os
REF=config["species"] + "." + config["genome"]
GENOME_VERSION = config["genome"]

rule assemble_transcripts:
    '''Rule adapted from ProteomeGenerator'''
    input:
        bam="{dir}/{sra}.sorted.bam" if check_sra() else "{dir}/{fq}.sorted.bam",
        gff=GFF3,
    output: "{dir}/{sra}.sorted.gtf" if check_sra() else "{dir}/{fq}.sorted.gtf"
    threads: 6
    log: "{dir}/{sra}.sorted.gtf.log" if check_sra() else "{dir}/{fq}.sorted.gtf"
    shell:
        "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -c 2.5 -m 300 -f .01 2> {log}" # strandedness: --fr for forwared or --rf for reverse

rule merge_transcripts:
    '''Rule adapted from ProteomeGenerator'''
    input:
        custom_gtfs=expand("{{dir}}/{sra}.sorted.gtf", sra=config["sra"]) if check_sra() is True else expand("{{dir}}/{fq}.sorted.gtf", fq=config["fq"]),
        gff=GFF3,
    output: "{dir}/combined.gtf"
    threads: 12
    benchmark: "{dir}/combined.gtf.benchmark"
    log: "{dir}/combined.gtf.log"
    shell:
        "stringtie --merge -o {output} -c 2.5 -m 300 -T 1 -f .01 -p {threads} -i {input.custom_gtfs} 2> {log}"

rule convert2ucsc:
    '''Rule adapted from ProteomeGenerator'''
    input: "{dir}/combined.gtf"
    output: "{dir}/combined_ucsc.gtf"
    shell: "python scripts/convert_ensembl2ucsc.py {input} {output}"

rule gtf_file_to_cDNA_seqs:
    '''Rule adapted from ProteomeGenerator'''
    input:
        fa=FA,
        gtf="{dir}/combined.gtf"
    output:
        fasta="{dir}/combined.transcripts.fasta",
        gtf="{dir}/combined.transcripts.gtf"
    benchmark: "{dir}/combined.gtf_file_to_cDNA_seqs.benchmark"
    log: "{dir}/combined.gtf_file_to_cDNA_seqs.log"
    threads: 1
    shell:
        "(gffread {input.gtf} -T -o {output.gtf} --no-pseudo --force-exons -M -Q && "
        "gffread -w {output.fasta} -g {input.fa} {output.gtf}) 2> {log}"

rule makeblastdb:
    '''Rule adapted from ProteomeGenerator'''
    input: UNIPROTFASTA
    output: [f"{UNIPROTFASTA}.pin", f"{UNIPROTFASTA}.phr", f"{UNIPROTFASTA}.psq"]
    benchmark: os.path.join(os.path.dirname(UNIPROTFASTA), os.path.basename(UNIPROTFASTA) + "makeblastdb.benchmark")
    log: os.path.join(os.path.dirname(UNIPROTFASTA), os.path.basename(UNIPROTFASTA) + "makeblastdb.log")
    threads: 1
    shell: "makeblastdb -in {input} -dbtype prot 2> {log}"

rule blastp:
    '''Rule adapted from ProteomeGenerator'''
    input:
        fasta=UNIPROTFASTA,
        pep="{dir}/longest_orfs.pep",
        blastdb=[UNIPROTFASTA+'.pin', UNIPROTFASTA+'.phr', UNIPROTFASTA+'.psq']
    output: "{dir}/combined.blastp.outfmt6"
    benchmark: "{dir}/combined.blastp.benchmark"
    log: "{dir}/combined.blastp.log"
    threads: 23
    shell: "blastp \
        -num_threads {threads} \
        -query {input.pep}  \
        -db {input.fasta}  \
        -max_target_seqs 1 \
        -outfmt 6 \
        -evalue 1e-5 \
        > {output} 2> {log}"

rule LongOrfs:
    '''Rule adapted from ProteomeGenerator'''
    input: "{dir}/combined.transcripts.fasta"
    output:
        "{dir}/longest_orfs.pep",
        # temp(directory("{dir}.__checkpoints_longorfs"))
    benchmark: "{dir}/combined.LongOrfs.benchmark"
    log: "{dir}/combined.LongOrfs.log"
    threads: 1
    shell: "TransDecoder.LongOrfs -O {wildcards.dir} -t {input} -m 100 2> {log}" # -S for strand-specific

rule Predict:
    '''Rule adapted from ProteomeGenerator'''
    input:
        orfs="{dir}/longest_orfs.pep",
        fasta="{dir}/combined.transcripts.fasta",
        blastp="{dir}/combined.blastp.outfmt6"
    output:
        "{dir}/combined.transcripts.fasta.transdecoder.pep",
        temp(directory("{dir}/..__checkpoints")),
        gff3="{dir}/combined.transcripts.fasta.transdecoder.gff3"
    benchmark: "{dir}/combined.Predict.benchmark"
    log: "{dir}/combined.Predict.log"
    threads: 1
    shell:
        "cd {wildcards.dir} && TransDecoder.Predict -O . -t ../{input.fasta} "
        "--single_best_only --retain_blastp_hits ../{input.blastp} 2> ../{log}"

rule gtf_to_alignment_gff3:
    '''Rule adapted from ProteomeGenerator'''
    input: "{dir}/combined.transcripts.gtf"
    output: "{dir}/combined.transcripts.gff3"
    benchmark: "{dir}/combined.gtf_to_alignment_gff3.benchmark"
    log: "{dir}/combined.gtf_to_alignment_gff3.log"
    threads: 1
    shell: "gtf_to_alignment_gff3.pl {input} > {output} 2> {log}"

rule cdna_alignment_orf_to_genome_orf:
    '''Rule adapted from ProteomeGenerator'''
    input:
        gff3="{dir}/combined.transcripts.gff3",
        fasta_td="{dir}/combined.transcripts.fasta",
        gff3_td="{dir}/combined.transcripts.fasta.transdecoder.gff3"
    output: "{dir}/combined.transcripts.genome.gff3"
    benchmark: "{dir}/combined.cdna_alignment_orf_to_genome_orf.benchmark"
    log: "{dir}/combined.cdna_alignment_orf_to_genome_orf.log"
    threads: 1
    shell: "cdna_alignment_orf_to_genome_orf.pl {input.gff3_td} {input.gff3} {input.fasta_td} > {output} 2> {log}"

rule gff3_file_to_bed:
    '''Rule adapted from ProteomeGenerator'''
    input: "{dir}/combined.transcripts.genome.gff3"
    output: "{dir}/combined.proteome.bed"
    benchmark: "{dir}/combined.gff3_file_to_bed.benchmark"
    log: "{dir}/combined.gff3_file_to_bed.log"
    threads: 1
    shell:
        "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_bed.pl /dev/stdin | tail -n +2 > {output} 2> {log}"

rule gff3_file_to_proteins:
    '''Rule adapted from ProteomeGenerator'''
    input: "{dir}/combined.transcripts.genome.gff3"
    output: "{dir}/combined.proteome.fasta"
    benchmark: "{dir}/combined.gff3_file_to_proteins.benchmark"
    log: "{dir}/combined.gff3_file_to_proteins.log"
    threads: 1
    shell: "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_proteins.pl --gff3 /dev/stdin --fasta {FA} | egrep -o '^[^*]+' > {output} 2> {log}"

rule reorderFASTA:
    '''Rule adapted from ProteomeGenerator
    Had to do this in Ubuntu 18.06: sudo ln -s /lib/x86_64-linux-gnu/libreadline.so.7 /lib/x86_64-linux-gnu/libreadline.so.6 '''
    input: "{dir}/combined.proteome.fasta"
    output: "{dir}/combined.proteome.unique.fasta"
    benchmark: "{dir}/combined.reorderFASTA.benchmark"
    log: "{dir}/combined.reorderFASTA.log"
    threads: 1
    script: "scripts/reorderFASTA.R"

rule remove_exon_and_utr_information:
    '''
    There is an issue in SnpEff where it assumes that CDS and exons end at the
    same position. When they do not, and when frame != 0 on the reverse strand,
    it breaks when trying to correct the frame. However, it works just fine if
    we remove the exons and just keep the CDS from TransDecoder.

    Minimal example for future reference:
    17	transdecoder	gene	77321555	77323854	.	-	.	ID=MSTRG.16910;Name=ORF%20type%3Acomplete%20len%3A149%20%28-%29%2Cscore%3D29.03
    17	transdecoder	mRNA	77321604	77323854	.	-	.	ID=MSTRG.16910.2.p3;Parent=MSTRG.16910;Name=ORF%20type%3Acomplete%20len%3A149%20%28-%29%2Cscore%3D29.03
    17	transdecoder	five_prime_UTR	77321604	77322851	.	-	.	ID=MSTRG.16910.2.p3.utr5p1;Parent=MSTRG.16910.2.p3
    17	transdecoder	exon	77321604	77322918	.	-	.	ID=MSTRG.16910.2.p3.exon1;Parent=MSTRG.16910.2.p3
    17	transdecoder	CDS	77322852	77322918	.	-	0	ID=cds.MSTRG.16910.2.p3;Parent=MSTRG.16910.2.p3
    17	transdecoder	exon	77323262	77323854	.	-	.	ID=MSTRG.16910.2.p3.exon2;Parent=MSTRG.16910.2.p3
    17	transdecoder	CDS	77323262	77323509	.	-	2	ID=cds.MSTRG.16910.2.p3;Parent=MSTRG.16910.2.p3
    17	transdecoder	three_prime_UTR	77323510	77323854	.	-	.	ID=MSTRG.16910.2.p3.utr3p1;Parent=MSTRG.16910.2.p3
    '''
    input: "{dir}/combined.transcripts.genome.gff3"
    output: "{dir}/combined.transcripts.genome.cds.gff3"
    log: "{dir}/combined.transcripts.genome.cds.log"
    shell: "python scripts/simplify_gff3.py {input} > {output} 2> {log}"

rule copy_gff3_to_snpeff:
    input: expand("{dir}/combined.transcripts.genome.cds.gff3", dir=config["analysisDirectory"])
    output: "SnpEff/data/combined.transcripts.genome.gff3/genes.gff"
    shell: "cp {input} {output}"

rule generate_snpeff_database:
    input:
        jar="SnpEff/snpEff.jar",
        gff3="SnpEff/data/combined.transcripts.genome.gff3/genes.gff",
        pfa="data/ensembl/" + REF + ".pep.all.fa",
        gfa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa"
    output:
        pfa="SnpEff/data/combined.transcripts.genome.gff3/protein.fa",
        gfa="SnpEff/data/genomes/combined.transcripts.genome.gff3.fa",
        done="SnpEff/data/combined.transcripts.genome.gff3/done.txt"
    params: ref="combined.transcripts.genome.gff3"
    resources: mem_mb=16000
    benchmark: "SnpEff/data/combined.transcripts.genome.gff3/snpeffdatabase.benchmark"
    log: "SnpEff/data/combined.transcripts.genome.gff3/snpeffdatabase.log"
    shell:
        "cp {input.pfa} {output.pfa} && "
        "cp {input.gfa} {output.gfa} && "
        "echo \"\n# {params.ref}\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.genome : Human genome {GENOME_VERSION} using RefSeq transcripts\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.M.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.MT.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "(java -Xmx{resources.mem_mb}M -jar {input.jar} build -gff3 -v {params.ref}) &> {log} && touch {output.done}"
