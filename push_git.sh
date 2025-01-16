#!/bin/bash

cur_dir=$(pwd)
cd /home/davidwang/hackweek2025/hubDirectory
git add .
git commit -m "."
git push https://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:@github.com/davidwang758/trackhub/ main
cd ${cur_dir}
