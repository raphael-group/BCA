#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob

my_path = "data/dynverse_data_v2/"
pca_files = glob.glob(my_path + '**/*.expression.outdir/*pca*.mtx', recursive=True)
bca_files = glob.glob(my_path + '**/*.expression.outdir/*bca*.mtx', recursive=True)
lda_files = glob.glob(my_path + '**/*.expression.outdir/*lda*.mtx', recursive=True)

# Rscript run_perf_structdr.R data/dynverse_data_v2/real/gold/aging-hsc-old_kowalczyk.expression.outidr/dim3_bca.niter30_ndim6_bw0.0_adaptive_bw10_automatic_bw2.0_relaxation0_ridge1_methodLocInv_stepsize0.5_maxangle90.0_k50.g_mst.mm.mtx

with open("3_commands_lda_structdr_perf.sh", 'w') as f:
    # for (i, file) in enumerate(pca_files):
    #     f.write(f"echo Running PCA-struct-perf on {file.split('/')[-1]}: {i+1}/339" + '\n')
    #     f.write(f"Rscript run_perf_structdr.R {file}" + '\n')


    # for (i, file) in enumerate(bca_files):
    #     f.write(f"echo Running BCA-struct on {file.split('/')[-1]}: {i+1}/339" + '\n')
    #     f.write(f"Rscript run_perf_structdr.R {file}" + '\n')


    for (i, file) in enumerate(lda_files):
        f.write(f"echo Running LDA-struct on {file.split('/')[-1]}: {i+1}/339" + '\n')
        f.write(f"Rscript run_perf_structdr.R {file}" + '\n')
