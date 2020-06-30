#!/bin/sh
cp /app/configs/config.yaml /app
read -r threads < /app/configs/threads.txt
snakemake -j $threads --restart-times 2
