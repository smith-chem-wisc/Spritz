# Spritz
Software for RNA-Seq analysis on Windows, including creating sample-specific proteoform databases from genomic data
[![Build status](https://ci.appveyor.com/api/projects/status/p54yrm6iixqm8jsf?svg=true)](https://ci.appveyor.com/project/acesnik/spritz)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c4187aae0b1f43c79e3e6379f77a408a)](https://www.codacy.com/app/acesnik/Spritz?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=smith-chem-wisc/Spritz&amp;utm_campaign=Badge_Grade)

Spritz uses the Windows Subsystem for Linux (WSL) to install and run commandline tools for Next-Generation Sequencing (NGS) analysis.

Spritz can be downloaded [here](https://github.com/smith-chem-wisc/Spritz/releases).

## Running Spritz with GUI

## Running Spritz without GUI

1. Launch Ubuntu and install snakemake: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html

2. Define config file as follows:

    #### Define SRA accession
    `sra: ["SRX159823"]`

    - Separate multiple SRAs with commas
    - Leave `sra: []` if using input fastqs.

    #### Define path to analysis directory
    `analysisDirectory: [/mnt/f/analysis]`

    - Follow this path format: /mnt/`{drive}`/`{folder}`

    #### Define species
    `species: "Homo_sapiens"`

    #### Define genome version
    `genome: "GRCh38"`

    #### Define Ensembl release
    `release: "81"`
	
	#### Define organism name
    `organism: "human"`
	
	- Based on Uniprot's naming convention
	- You can refer to genomes.csv

    #### DO NOT CHANGE: SnpEff version
    `snpeff: "86"`

    #### Define input fastq files paired
    `fq: [SRX159823]`

    - Fastq files should have *_1.fastq & *_2.fastq extensions and follow this naming convention: `{NAME}`_1.fastq and `{NAME}`_2.fastq
    - Leave `fq: []` if using SRA.

3. Run workflow:

    `snakemake -j 24`
