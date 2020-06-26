# Spritz
Software for RNA-Seq analysis on Windows, including creating sample-specific proteoform databases from genomic data

[![Build status](https://ci.appveyor.com/api/projects/status/p54yrm6iixqm8jsf?svg=true)](https://ci.appveyor.com/project/acesnik/spritz)
[![Release](https://img.shields.io/github/v/release/smith-chem-wisc/Spritz)](https://github.com/smith-chem-wisc/Spritz/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/smith-chem-wisc/Spritz/total.svg)](https://github.com/smith-chem-wisc/Spritz/releases/)
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

## Test it out! Try constructing the database for U2OS from the paper.

1. Add SRR629563 to the SRA list.

2. Create the Spritz workflow. Select "release-82" and "homo_sapiens."

3. Run Spritz!

Monitor progress by opening the `workflow_*.txt` file in [Notepad++](https://notepad-plus-plus.org/downloads/), which will ask to reload when updates have been made, i.e. when a new step has started or finished. The final database named `combined.spritz.snpeff.protein.withmods.xml.gz` can be used to search MS/MS with [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) to find variant peptides and proteoforms, possibly with modifications. We recommend performing 1) Calibration, 2) Global PTM Discovery (G-PTM-D), and 3) Search tasks to get the best results.

![image](https://user-images.githubusercontent.com/16342951/85874687-a76be700-b798-11ea-9bff-9f68646b03de.png)
