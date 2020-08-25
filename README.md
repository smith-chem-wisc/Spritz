# Spritz
Software for RNA-Seq analysis on Windows, including creating sample-specific proteoform databases from genomic data

[![Build status](https://ci.appveyor.com/api/projects/status/p54yrm6iixqm8jsf?svg=true)](https://ci.appveyor.com/project/acesnik/spritz)
[![Release](https://img.shields.io/github/v/release/smith-chem-wisc/Spritz)](https://github.com/smith-chem-wisc/Spritz/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/smith-chem-wisc/Spritz/total.svg)](https://github.com/smith-chem-wisc/Spritz/releases/)
[![Docker Pulls](https://img.shields.io/docker/pulls/smithlab/spritz)](https://hub.docker.com/r/smithlab/spritz/tags?page=1&ordering=last_updated)
[![Follow us on Twitter](https://img.shields.io/twitter/follow/smith_chem_wisc?label=Twitter&style=social)](https://twitter.com/smith_chem_wisc)

Spritz uses snakemake and Docker to install and run commandline tools for Next-Generation Sequencing (NGS) analysis.

Spritz can be downloaded [here](https://github.com/smith-chem-wisc/Spritz/releases).

![image](https://user-images.githubusercontent.com/16342951/84078314-55585280-a99e-11ea-9096-bebfcbb06bef.png)

## Running Spritz with GUI

1. Install [Docker Desktop for Windows](https://hub.docker.com/editions/community/docker-ce-desktop-windows).

2. Under Docker Settings, enable C and other necessary drives and allocate sufficient resources (recommended 16GB).

    ![settings](https://user-images.githubusercontent.com/42819128/70090841-8a937a80-15e0-11ea-9742-ca959a89deca.png)

3. Launch Spritz.

    ![howto](https://user-images.githubusercontent.com/42819128/70091146-2624eb00-15e1-11ea-9230-bfd118aa03d9.png)

    Step 1: Input SRA accessions OR FASTQ files.

    Step 2: Create and customize your Spritz workflow.
    
   ![workflow](https://user-images.githubusercontent.com/42819128/70091992-e65f0300-15e2-11ea-9e0f-7bb4262afefa.png)
   
    Step 3: Run Spritz!

## System Requirements

* Environment:
  * .NET Core 3.1:
     * Windows: https://dotnet.microsoft.com/download/dotnet-core/thank-you/runtime-desktop-3.1.3-windows-x64-installer
* Note that the installer (MetaMorpheusInstaller.msi) only works on Windows. The command-line version of MetaMorpheus supports any operating system that supports .NET Core (Windows, MacOS, Linux)
* 16 GB RAM recommended

## Test it out! Try constructing the database for U2OS from the paper.

1. Add SRR629563 to the SRA list.

2. Create the Spritz workflow. Select "release-82" and "homo_sapiens."

3. Run Spritz!

Monitor progress in the Information textbox. The final database named `final/combined.spritz.snpeff.protein.withmods.xml.gz` can be used to search MS/MS with [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) to find variant peptides and proteoforms, possibly with modifications. We recommend performing 1) Calibration, 2) Global PTM Discovery (G-PTM-D), and 3) Search tasks to get the best results.

![image](https://user-images.githubusercontent.com/16342951/85874687-a76be700-b798-11ea-9bff-9f68646b03de.png)

The final database named `final/combined.spritz.snpeff.protein.fasta` is generated to contain variant protein sequences, and it may be used in other search software, such as Proteome Discoverer, ProSight, and MASH Explorer.

The final database named `final/combined.spritz.snpeff.protein.withdecoys.fasta` is ready for use in MSFragger. It is generated to contain variant protein sequences with decoy protein sequences appended.
