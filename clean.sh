#!/bin/bash

error=true
until [ "$error" = false ]; do 
	(gatk IndexFeatureFile -F data/ensembl/Homo_sapiens.clean.vcf) &> vcf_index.txt	
	line=$(grep 'line number' vcf_index.txt)
	if [ -z "$line" ]
	then 
		error=false
	else
		num=$(echo $line | tr -dc '0-9')
		sed -i ${num}d data/ensembl/Homo_sapiens.clean.vcf
	fi
	rm vcf_index.txt
done