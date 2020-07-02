# script to get download ensembl references
# assumes all paths exist (bash script in gui)

import subprocess
import sys
import os

def check_url_exists(url):
    return subprocess.check_output(['./validate.sh', url]) == b'true\n'

speciesOrig = sys.argv[1].split('.')[0]
species = speciesOrig.lower()
path = '.'.join(sys.argv[1].split('.')[:-1]) # Mus_musculus.GRCm38.86
release = sys.argv[1].split('.')[-1]
fileToGet = sys.argv[2] # vcf, not, or dry
protocol = "http" #ftp or http; http seems to work in more places

primary = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/dna/{path}.dna.primary_assembly.fa.gz"
toplevel = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/dna/{path}.dna.toplevel.fa.gz"

vcf1 = f"{protocol}://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{speciesOrig}.vcf.gz"
vcf2 = f"{protocol}://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{species}.vcf.gz" # edit, incosistent naming convention /vcf/Mus_musculus.vcf.gz or /vcf/mus_musculus.vcf.gz
gff = f"{protocol}://ftp.ensembl.org/pub/release-{release}/gff3/{species}/{path}.{release}.gff3.gz"
pep = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/pep/{path}.pep.all.fa.gz"

if fileToGet == "vcf":
    if check_url_exists(vcf1):
        subprocess.check_call(["wget", "-P", "data/ensembl/", vcf1])
    else:
        subprocess.check_call(["wget", "-P", "data/ensembl/", vcf2])
        os.rename(f"data/ensembl/{species}.vcf.gz", f"data/ensembl/{speciesOrig}.vcf.gz")

elif fileToGet == "not":
    if check_url_exists(primary):
        subprocess.check_call(["wget", "-P", "data/ensembl/", primary])
    else:
        subprocess.check_call(["wget", "-P", "data/ensembl/", toplevel])
    subprocess.check_call(["wget", "-P", "data/ensembl/", gff])
    subprocess.check_call(["wget", "-P", "data/ensembl/", pep])

else:
    # perform a dry run of the download_sras
    sys.stdout.write(f"{path}.{release}\tprimary\t{check_url_exists(primary)}\t{primary}\n")
    sys.stdout.write(f"{path}.{release}\ttoplevel\t{check_url_exists(toplevel)}\t{toplevel}\n")
    sys.stdout.write(f"{path}.{release}\tvcf1\t{check_url_exists(vcf1)}\t{vcf1}\n")
    sys.stdout.write(f"{path}.{release}\tvcf2\t{check_url_exists(vcf2)}\t{vcf2}\n")
    sys.stdout.write(f"{path}.{release}\tgff\t{check_url_exists(gff)}\t{gff}\n")
    sys.stdout.write(f"{path}.{release}\tpep\t{check_url_exists(pep)}\t{pep}\n")
