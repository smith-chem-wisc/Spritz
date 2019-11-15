import sys
import numpy as np
import pandas as pd
import os
import glob

USAGE = "python make_rsem_dataframe.py <column_index> <gff> <outCountsFile> <outTpmsFile> <namesOut>"
if len(sys.argv) != 5: print(USAGE); exit();
gff, outcounts, outtpms, outn = sys.argv[1:]

print(f"getting gene names from {gff}")
with open(gff) as gffhandle:
    line_ct = sum([1 for line in gffhandle])
gene_id_to_name, gene_id_to_biotype, gene_id_to_description = {}, {}, {}
with open(gff) as gffhandle:
    for i, line in enumerate(gffhandle):
        if i % 50000 == 0: print(f"Processed {i} lines out of {line_ct} from {gff}.")
        if not line.startswith('#'):
            line = line.strip().rstrip(';').split('\t')
            if len(line) < 9: continue
            attribute_list = line[8].replace('; ', ';').split(';')
            attributes = {}
            for item in attribute_list:
                item = item.split('=')
                attributes[item[0]] = item[1]
            if "ID" in attributes and "Name" in attributes and attributes["ID"].startswith("gene:"):
                gene_id_to_name[attributes["ID"][len("gene:"):]] = attributes["Name"]
            else: continue
            gene_id_to_biotype[attributes["ID"][len("gene:"):]] = attributes["biotype"] if "ID" in attributes and "biotype" in attributes else ""
            gene_id_to_description[attributes["ID"][len("gene:"):]] = attributes["description"] if "ID" in attributes and "description" in attributes else ""


print("globbing output/*.genes.results")
files = glob.glob("output/*.genes.results")
line_ct = sum(1 for line in open(files[0]))

ids = [line.split('\t')[0].strip("gene:") for line in open(files[0])]
ids[0] = "gene_id"
names = [gene_id_to_name[id] for id in ids[1:]]
biotypes = [gene_id_to_biotype[id] for id in ids[1:]]
description = [gene_id_to_description[id] for id in ids[1:]]
counts = [np.array(ids)]
tpms = [np.array(ids)]
for file in files:
    counts.append(np.array([line.split('\t')[4] for line in open(file)]))
    tpms.append(np.array([line.split('\t')[5] for line in open(file)]))
    print("reading "+ os.path.basename(file).split(".")[0])

print(f"Saving to {outcounts} ...")
dataframe = np.row_stack(counts)
dataframe[1:,0] = [os.path.basename(file).split(".")[0]]
pddf = pd.DataFrame(dataframe[1:,1:], index=dataframe[1:,0], columns=dataframe[0,1:]).sort_index()
pddf.filter(regex="ENSG").to_csv(outcounts)
pddf.filter(regex="ERCC").to_csv(outcounts + ".ercc.csv")

print(f"Saving to {outtpms} ...")
dataframe = np.row_stack(tpms)
dataframe[1:,0] = [os.path.basename(file).split(".")[0] for file in files]
pddf = pd.DataFrame(dataframe[1:,1:], index=dataframe[1:,0], columns=dataframe[0,1:]).sort_index()
pddf.filter(regex="ENSG").to_csv(outtpms)
pddf.filter(regex="ERCC").to_csv(outtpms + ".ercc.csv")

print(f"Saving to {outn} ...")
np.savetxt(outn, np.column_stack((ids[1:], names, biotypes, description)), delimiter=",", fmt="%s")

print("Done.")
