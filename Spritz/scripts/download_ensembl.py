# script to get download ensembl references
# assumes all paths exist (bash script in gui)

import yaml
import subprocess
import sys
import os

with open("config.yaml", 'r') as stream:
    data = yaml.safe_load(stream)

species = data["species"].lower()
release = data["release"]

path = sys.argv[1] # Mus_musculus.GRCm38
fileToGet = sys.argv[2] # vcf or not
protocol = "http" #ftp or http; http seems to work in more places

primary = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/dna/{path}.dna.primary_assembly.fa.gz"
toplevel = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/dna/{path}.dna.toplevel.fa.gz"

if fileToGet != "vcf":
    # download gfa
    if subprocess.check_output(['./validate.sh', primary]) == b'true\n':
        subprocess.check_call(["wget", "-P", "data/ensembl/", primary])
    else:
        subprocess.check_call(["wget", "-P", "data/ensembl/", toplevel])

    # download gff
    gff = f"{protocol}://ftp.ensembl.org/pub/release-{release}/gff3/{species}/{path}.{release}.gff3.gz"
    subprocess.check_call(["wget", "-P", "data/ensembl/", gff])

    # download pep
    pep = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/pep/{path}.pep.all.fa.gz"
    subprocess.check_call(["wget", "-P", "data/ensembl/", pep])

else:
    # download vcf
    vcf1 = f"{protocol}://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{data['species']}.vcf.gz"
    vcf2 = f"{protocol}://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{species}.vcf.gz" # edit, incosistent naming convention /vcf/Mus_musculus.vcf.gz or /vcf/mus_musculus.vcf.gz
    if subprocess.check_output(['./validate.sh', vcf1]) == b'true\n':
        subprocess.check_call(["wget", "-P", "data/ensembl/", vcf1])
    else:
        subprocess.check_call(["wget", "-P", "data/ensembl/", vcf2])
        os.rename(f"data/ensembl/{species}.vcf.gz", f"data/ensembl/{data['species']}.vcf.gz")
