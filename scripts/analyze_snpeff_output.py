# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 18:32:13 2019

@author: Anthony
"""

import argparse
from lxml import etree
import numpy as np
import pandas as pd
from os import path

#%%

argparser = argparse.ArgumentParser()
argparser.add_argument("--snpeff_html", "-s", help="SnpEff HTML to analyze", required=True)
args = argparser.parse_args()
snpeff_html = args.snpeff_html

#%%
# Idea: take
snpeff_html = "E:\\source\\repos\\SpritzSnake\\data\\combined.spritz.snpeff.html"
folder = path.dirname(snpeff_html)
prefix = path.basename(snpeff_html)[:-5]
htmlparser = etree.HTMLParser()
tree = etree.parse(snpeff_html, htmlparser)

#%%
body = tree.getroot()[1]
summary, changesByType, effects, baseChanges, alleleFreq,codonChanges,aaChanges = [],[],[],[],[],[],[],[]
for item in body.getiterator("a"):
    if item.get("name") == "summary": summary = item
    elif item.get("name") == "changesByType": changesByType = item
    elif item.get("name") == "effects": effects = item
    elif item.get("name") == "baseChanges": baseChanges = item
    elif item.get("name") == "alleleFreq": alleleFreq = item
    elif item.get("name") == "codonChanges": codonChanges = item
    elif item.get("name") == "aaChanges": aaChanges = item
    
#%% Summary information
genome, numvariants, genomelength, variantRate = "","","",""
for item in summary.getiterator("b"):
    label = item.text.strip()
    if label == "Genome": 
        genome = item.getparent().getnext().text.strip()
    elif label == "Number of variants (before filter)": 
        numvariants = item.getparent().getnext().text.strip()
    elif label == "Genome effective length": 
        variantRate = item.getparent().getnext().text.strip()
    elif label == "Variant rate": 
        variantRate = item.getparent().getnext().text.strip()

print(f"Genome: {genome}")
print(f"Number of Variants: {numvariants}")
print(f"Genome Length: {genomelength}")
print(f"Variant Rate: {variantRate}")

#%% Changes by type
vartypes = []
varcount = []
for item in changesByType.getiterator("b"):
    if item.getparent().getparent().getparent().tag == "tbody":
        label = item.text.strip()
        num = item.getparent().getnext().text.strip()
        vartypes.append(label)
        varcount.append(num)
        
changesByTypeDf = pd.DataFrame({"Type" : vartypes, "Total": varcount})
print(changesByTypeDf.to_string(index=False))
changesByTypeDf.to_csv(path.join(folder, f"{prefix}.changesByType.csv"), index=False)

#%% Effects
efftypes = []
effcount = []
for item in effects.getiterator("b"):
    if item.getparent().tag == "td":
        label = item.text.strip()
        num = item.getparent().getnext().getnext().text.strip()
        efftypes.append(label)
        effcount.append(int(num.replace(",","")))
        
efftypes, effcount = np.array(efftypes), np.array(effcount)
effloc_idx = [all([ee.isupper() for ee in eff]) for eff in efftypes]
changesByTypeDf = pd.DataFrame({"Type" : efftypes, "Total": effcount})
print(changesByTypeDf.to_string(index=False))
changesByTypeDf.to_csv(path.join(folder, f"{prefix}.effects.csv"), index=False)

# could print synonymous_variant, missense_variant, frameshift_variant, stop_gained, 
# conservative_inframe_insertion + disruptive_inframe_insertion,
# conservative_inframe_deletion + disruptive_inframe_deletion, stop_loss, since those are the ones used
# in the proteogenomics

#%% base changes
#%% allele frequencies
#%% codon changes
#%% amino acid changes

