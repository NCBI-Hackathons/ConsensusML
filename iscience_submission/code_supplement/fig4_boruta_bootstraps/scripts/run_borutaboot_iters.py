#!/usr/bin/env python3

# run iterations of borutaboot_1ter

import subprocess

iterindices = [ri for ri in range(1,1001)]

for ii, ri in enumerate(iterindices, 0):
    seed = ii
    print("beginning iter "+str(ii))
    subprocess.call(['Rscript','borutaboot_iter_nrank.R',"nrank_"+str(ri),str(seed)], shell=False)
    print("finished iter "+str(ii))