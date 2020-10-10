# Spritz
Software for RNA-Seq analysis on Windows, including creating sample-specific proteoform databases from genomic data

[![Build status](https://ci.appveyor.com/api/projects/status/p54yrm6iixqm8jsf?svg=true)](https://ci.appveyor.com/project/acesnik/spritz)
[![Release](https://img.shields.io/github/v/release/smith-chem-wisc/Spritz)](https://github.com/smith-chem-wisc/Spritz/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/smith-chem-wisc/Spritz/total.svg)](https://github.com/smith-chem-wisc/Spritz/releases/)
[![Docker Pulls](https://img.shields.io/docker/pulls/smithlab/spritz)](https://hub.docker.com/r/smithlab/spritz/tags?page=1&ordering=last_updated)
[![Follow us on Twitter](https://img.shields.io/twitter/follow/smith_chem_wisc?label=Twitter&style=social)](https://twitter.com/smith_chem_wisc)

Spritz uses snakemake and Docker to install and run commandline tools for Next-Generation Sequencing (NGS) analysis.

Spritz can be downloaded [here](https://github.com/smith-chem-wisc/Spritz/releases).

![image](https://user-images.githubusercontent.com/16342951/93618988-a3b5be00-f99d-11ea-8be4-063395e24ce1.png)

## Running Spritz with GUI

1. Install [Docker Desktop for Windows](https://hub.docker.com/editions/community/docker-ce-desktop-windows).

2. Allocate resources to Docker. 

    This is the preferred method; it requires Windows 10 version 2004, but it is more robust. Proceed to the alternate method below if you do not have Windows version 2004 or greater installed. 

    This utilizes a new Windows feature, [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install-win10), that provides a native Linux subsystem that [Docker can utilize](https://docs.docker.com/docker-for-windows/wsl/).

    Steps displayed below:
    * Check that your Windows version is > 2004 in "About your PC"
    * Install Ubuntu 20.04 from Microsoft store
    * Enable the WSL feature
    * Restart your computer
    * Start up Ubuntu 20.04 and complete setup
    * Ensure that WSL2 is enabled in Docker Desktop.

    ![wsl2setup](https://user-images.githubusercontent.com/16342951/95661576-03633d00-0af6-11eb-90e5-2f12cf9f5c18.png)
    
    This is the alternate method, which is less robust but available on all Windows versions.

    This utilizes the virtual machine "Vmmem" that comes with Docker.

    First, ensure that the C:\ drive and any data drive you plan to use are enabled for use. 
    
    Then, allocate enough resources: 
    * Allow all CPUs to be utilized.
    * Provide at least 16 MB of RAM, but more is better, and allowing Docker to utilize all available RAM is most robust.
    * Provide at least 60 GB of disk space, but more is better, and allowing Docker to utilize all available disk space is most robust.

    ![vmmemSetup](https://user-images.githubusercontent.com/16342951/95661738-3c4fe180-0af7-11eb-8d61-a2c86172a22c.png)

4. Launch Spritz.

    ![howto](https://user-images.githubusercontent.com/42819128/70091146-2624eb00-15e1-11ea-9230-bfd118aa03d9.png)

    Step 1: Input SRA accessions OR FASTQ files.

    Step 2: Create and customize your Spritz workflow.

   ![workflow](https://user-images.githubusercontent.com/42819128/70091992-e65f0300-15e2-11ea-9e0f-7bb4262afefa.png)

    Step 3: Run Spritz!

## System Requirements

* Environment:
  * .NET Core 3.1:
     * Windows: https://dotnet.microsoft.com/download/dotnet-core/thank-you/runtime-desktop-3.1.3-windows-x64-installer
* 16 GB RAM recommended
* Note that the installer (Spritz.msi) only works on Windows. 
* Spritz will also work on the commandline within a Unix system (Linux, Mac, WSL on Windows). First, install [miniconda3](https://docs.conda.io/en/latest/miniconda.html), and then create a `conda` environment for spritz by running `conda env create --name spritz --file environment.yaml; conda activate spritz`. After adapting the `config.yaml` file manually, Spritz may be run using `snakemake -j {threads} --resources mem_mb={memory_megabytes}`, where `{threads}` and `{memory_megabytes}` are replaced with your specifications.

## Test it out! Try constructing the database for U2OS from the paper.

1. Add SRR629563 to the SRA list.

2. Create the Spritz workflow. Select "release-82" and "homo_sapiens."

3. Run Spritz!

Monitor progress in the Information textbox. The final database named `final/combined.spritz.snpeff.protein.withmods.xml.gz` can be used to search MS/MS with [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) to find variant peptides and proteoforms, possibly with modifications. We recommend performing 1) Calibration, 2) Global PTM Discovery (G-PTM-D), and 3) Search tasks to get the best results.

![image](https://user-images.githubusercontent.com/16342951/85874687-a76be700-b798-11ea-9bff-9f68646b03de.png)

The final database named `final/combined.spritz.snpeff.protein.fasta` is generated to contain variant protein sequences, and it may be used in other search software, such as Proteome Discoverer, ProSight, and MASH Explorer.

The final database named `final/combined.spritz.snpeff.protein.withdecoys.fasta` is ready for use in MSFragger. It is generated to contain variant protein sequences with decoy protein sequences appended.
