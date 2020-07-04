# script to get download ensembl references
# assumes all paths exist (bash script in gui)

import subprocess
import sys
import os
import requests

def check_url_exists(url):
    return requests.get(url).status_code == 200

speciesOrig = sys.argv[1].split('.')[0]
species = speciesOrig.lower()
path = '.'.join(sys.argv[1].split('.')[:-1]) # Mus_musculus.GRCm38.86
release = sys.argv[1].split('.')[-1]
fileToGet = sys.argv[2] # vcf, not, or dry
protocol = "http" #ftp or http; http seems to work in more places

primary = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/dna/{path}.dna.primary_assembly.fa.gz"
toplevel = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/dna/{path}.dna.toplevel.fa.gz"
gff = f"{protocol}://ftp.ensembl.org/pub/release-{release}/gff3/{species}/{path}.{release}.gff3.gz"
pep = f"{protocol}://ftp.ensembl.org/pub/release-{release}//fasta/{species}/pep/{path}.pep.all.fa.gz"

if fileToGet == "not":
    gfa = primary if check_url_exists(primary) else toplevel
    subprocess.check_call(["wget", "-P", "data/ensembl/", gfa])
    subprocess.check_call(["wget", "-P", "data/ensembl/", gff])
    subprocess.check_call(["wget", "-P", "data/ensembl/", pep])

else:
    # perform a dry run of the download
    vcf1 = f"http://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{speciesOrig}.vcf.gz"
    vcf2 = f"http://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{species}.vcf.gz"
    sys.stdout.write(f"{path}.{release}\tprimary\t{check_url_exists(primary)}\t{primary}\n")
    sys.stdout.write(f"{path}.{release}\ttoplevel\t{check_url_exists(toplevel)}\t{toplevel}\n")
    sys.stdout.write(f"{path}.{release}\tgff\t{check_url_exists(gff)}\t{gff}\n")
    sys.stdout.write(f"{path}.{release}\tpep\t{check_url_exists(pep)}\t{pep}\n")
    sys.stdout.write(f"{path}.{release}\tgff\t{check_url_exists(vcf1)}\t{vcf1}\n")
    sys.stdout.write(f"{path}.{release}\tpep\t{check_url_exists(vcf2)}\t{vcf2}\n")
