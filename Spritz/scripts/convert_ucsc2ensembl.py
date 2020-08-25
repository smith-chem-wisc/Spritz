ucsc=open("data/ensembl/common_all_20170710.vcf")

ucsc2ensembl={}
for line in open("ChromosomeMappings/GRCh38_UCSC2ensembl.txt"):
    linesplit=line.strip().split("\t")
    if len(linesplit) <= 1: continue
    ucsc2ensembl[linesplit[0]] = linesplit[1]

with open("data/ensembl/common_all_20170710.ensembl.vcf","w") as ensembl:
    for line in ucsc:
        if line.startswith("#"):
            ensembl.write(line)
        splitline = line.split("\t")
        if len(splitline) > 1 and splitline[0] in ucsc2ensembl:
            splitline[0] = ucsc2ensembl[splitline[0]]
            ensembl.write("\t".join(splitline))
