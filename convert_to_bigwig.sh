#!/bin/bash

files=($(ls /home/davidwang/hackweek2025/hubDirectory/${1}/data | grep ".wig"))
for f in "${files[@]}"; do
	root=$(basename ${f} .wig)
	/home/davidwang/hackweek2025/./wigToBigWig /home/davidwang/hackweek2025/hubDirectory/${1}/data/${root}.wig /home/davidwang/hackweek2025/hg38.chrom.sizes /home/davidwang/hackweek2025/hubDirectory/${1}/data/${root}.bw
done
