import yaml
import re

with open("config.yaml", 'r') as stream:
   data = yaml.safe_load(stream)

species = data["species"]
version = data["genome"]

vcf=open("./data/ensembl/" + species + ".vcf")

def invalid(str):
    """ Check whether sequence str contains ANY of the items in set. """
    set = ['A','C','T','G']
    return True not in [c in str for c in set]

with open("./data/ensembl/" + species + ".clean.vcf","w") as ensembl:
    for line in vcf:
        # header
        if line.startswith("#"):
            ensembl.write(line)
            continue

        # remove any lines with empty alleles
        splitline = line.split("\t")
        if '' in splitline[0:7]: continue

        # remove unparsable alleles
        if invalid(splitline[4]): continue

        # remove duplicate alleles
        duplicate = re.findall("([a-zA-z]\.)", splitline[4])
        if len(duplicate) > 1: continue

        ensembl.write(line)
