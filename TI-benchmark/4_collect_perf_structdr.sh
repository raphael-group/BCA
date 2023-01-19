# rm structdr.perf.txt;
rm bca.structdr.perf.txt; for i in `find ./ -name '*_bca*perf.prune.txt'`;do cat $i |tail -n +2 >> bca.structdr.perf.txt;done
for i in `find ./ -name '*_lda*perf.prune.txt'`;do cat $i |tail -n +2 >> lda.structdr.perf.txt;done
for i in `find ./ -name '*_pca*perf.prune.txt'`;do cat $i |tail -n +2 >> pca.structdr.perf.txt;done
