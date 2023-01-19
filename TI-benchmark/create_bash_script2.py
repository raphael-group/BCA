#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob

my_path = "data/dynverse_data_v2/"
pca_files = glob.glob(my_path + '**/*.expression.pca*', recursive=True)
bca_files = glob.glob(my_path + '**/*.expression.bca*', recursive=True)
lda_files = glob.glob(my_path + '**/*.expression.lda*', recursive=True)

# python3 quasildr/run_structdr.py --input ./data/dynverse_data_v2/real/gold/aging-hsc-old_kowalczyk.expression.bca3 --output ./data/dynverse_data_v2/real/gold/aging-hsc-old_kowalczyk.expression.outidr/dim3_bca --ndim 6 --n_jobs 1 --anno_file ./data/dynverse_data_v2/real/gold/aging-hsc-old_kowalczyk.anno --extract_backbone

with open("2_commands_our_structdr.sh", 'w') as f:
    # for (i, file) in enumerate(pca_files):
    #     # The following line will extract the num_dim from embedding files
    #     # it isn't super robust but for our use case it is fine, taken
    #     # from https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python
    #     num_dim = int(''.join(filter(str.isdigit, file.split('.')[-1])))

    #     f.write(f"echo Running PCA-struct on {file.split('/')[-1]}: {i+1}/339" + '\n')

    #     prefix = "python3 quasildr/run_structdr.py --input "
    #     base = f"{file} --output {'.'.join(file.split('.')[:-1])}.outdir/dim{num_dim}_pca "
    #     suffix = f"--ndim 6 --n_jobs 1 --anno_file {file.split('.')[0]}.anno --extract_backbone"
    #     f.write(prefix + base + suffix + '\n')


    for (i, file) in enumerate(bca_files):
        num_dim = int(''.join(filter(str.isdigit, file.split('.')[-1])))
        f.write(f"echo Running BCA-struct on {file.split('/')[-1]}: {i+1}/339" + '\n')

        prefix = "python3 quasildr/run_structdr.py --input "
        base = f"{file} --output {'.'.join(file.split('.')[:-1])}.outdir/dim{num_dim}_bca "
        suffix = f"--ndim 6 --n_jobs 1 --anno_file {file.split('.')[0]}.anno --extract_backbone"
        f.write(prefix + base + suffix + '\n')


    for (i, file) in enumerate(lda_files):
        num_dim = int(''.join(filter(str.isdigit, file.split('.')[-1])))
        f.write(f"echo Running LDA-struct on {file.split('/')[-1]}: {i+1}/339" + '\n')

        prefix = "python3 quasildr/run_structdr.py --input "
        base = f"{file} --output {'.'.join(file.split('.')[:-1])}.outdir/dim{num_dim}_lda "
        suffix = f"--ndim 6 --n_jobs 1 --anno_file {file.split('.')[0]}.anno --extract_backbone"
        f.write(prefix + base + suffix + '\n')
