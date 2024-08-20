#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import Levenshtein as lv

# %matplotlib inline

#pairs of R1 files
liblist1 = ["/home/jeannie/group/OvCa_Genomic_Data/newCloneCalling/allowlisting_scripts/H180_2_scRNAseq_S2_R1_celltag_reads_dist_3.txt"]
liblist2 = ["/home/jeannie/group/OvCa_Genomic_Data/newCloneCalling/allowlisting_scripts/H180_1_scRNAseq_S1_R1_celltag_reads_01_dist_3.txt"]
liblist3 = ["/home/jeannie/group/OvCa_Genomic_Data/newCloneCalling/allowlisting_scripts/H180_1_scRNAseq_S1_R1_celltag_reads_02_dist_3.txt"]
# split sample 1 into 2 halves, then add liblist3
lib_list = []
celltags = {}

for i, j, k in zip(liblist1,liblist2, liblist3):
    lib1 = pd.read_csv(i, sep="\t", header = None)
    lib2 = pd.read_csv(j, sep="\t", header = None)
    lib3 = pd.read_csv(k, sep="\t", header = None)
    lib1['sum'] = lib1[1].cumsum()
    lib2['sum'] = lib2[1].cumsum()
    lib3['sum'] = lib3[1].cumsum()
    th_curr = max(10, np.ceil(np.mean(np.array([np.percentile(lib1[1],90)/10, np.percentile(lib2[1],90)/10, np.percentile(lib3[1],90)/10]))))
    print(f'Read threshold used: {th_curr}')
    
    lib1_fil = lib1[lib1[1]>th_curr]
    lib2_fil = lib2[lib2[1]>th_curr]
    lib3_fil = lib3[lib3[1]>th_curr]
    # lib_curr = set(lib1_fil[0]).intersection(set(lib2_fil[0]))

    sample1 = set(lib2_fil[0]) | set(lib3_fil[0])
    lib_curr = set(lib1_fil[0]) & sample1
    #lib_curr = set(lib1_fil[0]) & set(lib2_fil[0]) & set(lib3_fil[0])
    lib_list.append(lib_curr)

final_lib = list(set().union(*lib_list))
pd.DataFrame(final_lib).to_csv("allowlist.csv", index = False, header=False)