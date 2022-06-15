#!/usr/bin/bash

echo "Starting Condor"
/home/egalcera/anaconda3/envs/tfgpu/bin/python JANE-gpu.py --no-confirm --config-file config2.txt
echo "Finishing Condor"
