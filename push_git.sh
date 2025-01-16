#!/bin/bash

cur_dir=$(pwd)
cd /home/davidwang/hackweek2025/hubDirectory
git add .
git commit -m "."
#have to remove the push because it contains my access token
cd ${cur_dir}
