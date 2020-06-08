# Spritz
Software for RNA-Seq analysis on Windows, including creating sample-specific proteoform databases from genomic data
[![Build status](https://ci.appveyor.com/api/projects/status/p54yrm6iixqm8jsf?svg=true)](https://ci.appveyor.com/project/acesnik/spritz)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c4187aae0b1f43c79e3e6379f77a408a)](https://www.codacy.com/app/acesnik/Spritz?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=smith-chem-wisc/Spritz&amp;utm_campaign=Badge_Grade)

Spritz uses the Windows Subsystem for Linux (WSL) to install and run commandline tools for Next-Generation Sequencing (NGS) analysis.

Spritz can be downloaded [here](https://github.com/smith-chem-wisc/Spritz/releases).

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
