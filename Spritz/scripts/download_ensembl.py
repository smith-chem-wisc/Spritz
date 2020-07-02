# script to get download ensembl references
# assumes all paths exist (bash script in gui)

import yaml
import subprocess
import sys
import os

def check_url_exists(url):
    return subprocess.check_output(['./validate.sh', url]) == b'true\n'

with open("config.yaml", 'r') as stream:
    data = yaml.safe_load(stream)

species = data["species"].lower()
release = data["release"]

path = sys.argv[1] # Mus_musculus.GRCm38
fileToGet = sys.argv[2] # vcf, not, or dry
protocol = "http" #ftp or http; http seems to work in more places

primary = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/dna/{path}.dna.primary_assembly.fa.gz"
toplevel = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/dna/{path}.dna.toplevel.fa.gz"

vcf1 = f"{protocol}://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{data['species']}.vcf.gz"
vcf2 = f"{protocol}://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{species}.vcf.gz" # edit, incosistent naming convention /vcf/Mus_musculus.vcf.gz or /vcf/mus_musculus.vcf.gz
gff = f"{protocol}://ftp.ensembl.org/pub/release-{release}/gff3/{species}/{path}.{release}.gff3.gz"
pep = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/pep/{path}.pep.all.fa.gz"

if fileToGet == "vcf":
    if check_url_exists(vcf1):
        subprocess.check_call(["wget", "-P", "data/ensembl/", vcf1])
    else:
        subprocess.check_call(["wget", "-P", "data/ensembl/", vcf2])
        os.rename(f"data/ensembl/{species}.vcf.gz", f"data/ensembl/{data['species']}.vcf.gz")

elif fileToGet == "not":
    if check_url_exists(primary):
        subprocess.check_call(["wget", "-P", "data/ensembl/", primary])
    else:
        subprocess.check_call(["wget", "-P", "data/ensembl/", toplevel])
    subprocess.check_call(["wget", "-P", "data/ensembl/", gff])
    subprocess.check_call(["wget", "-P", "data/ensembl/", pep])

else:
    # perform a dry run of the download_sras
    sys.stdout.write(f"{path}\tprimary\t{check_url_exists(primary)}\t{primary}\n")
    sys.stdout.write(f"{path}\ttoplevel\t{check_url_exists(toplevel)}\t{toplevel}\n")
    sys.stdout.write(f"{path}\tvcf1\t{check_url_exists(vcf1)}\t{vcf1}\n")
    sys.stdout.write(f"{path}\tvcf2\t{check_url_exists(vcf2)}\t{vcf2}\n")
    sys.stdout.write(f"{path}\tgff\t{check_url_exists(gff)}\t{gff}\n")
    sys.stdout.write(f"{path}\tpep\t{check_url_exists(pep)}\t{pep}\n")
