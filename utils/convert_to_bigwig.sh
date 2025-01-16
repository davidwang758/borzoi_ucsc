#!/bin/bash

files=($(ls ${2}${1}/data | grep ".wig"))
for f in "${files[@]}"; do
	root=$(basename ${f} .wig)
	${3}./wigToBigWig ${2}${1}/data/${root}.wig ${3}hg38.chrom.sizes ${2}${1}/data/${root}.bw
done
