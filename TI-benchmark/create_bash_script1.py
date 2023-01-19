#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob

my_path = "data/dynverse_data_v2/"
# files = glob.glob(my_path + '**/*.expression', recursive=True)
files = glob.glob(my_path + '**/*.expression.bca1', recursive=True)

with open("1_commands_our_dr.sh", 'w') as f:
    for (i, file) in enumerate(files):
        f.write(f"echo Running {file.split('/')[-1]}: {i+1}/339" + '\n')
        f.write(f"python3 run_pca_lda_bca.py {file}" + '\n')
