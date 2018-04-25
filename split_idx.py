#usage: generate the 6 snp-index files geno1-geno6, prepare pheno, parents.txt, sample_idx with a train/test split.
#Then, do:
#
#python split_idx.py parents.txt sample_idx > sample_idx_with_validation
#
#Now do:
#
#python brr_val.py geno1 pheno parents.txt sample_idx_with_validation 2
#python brr_val.py geno2 pheno parents.txt sample_idx_with_validation 2
#...

import csv
import sys
import random
import scipy as SP

n_s=6210
parents = SP.array(list(csv.reader(open(sys.argv[1],'rb'),
                                   delimiter='\t'))).astype(int)
split_idx=list(csv.reader(open(sys.argv[2],'rb'),delimiter='\t'))
train=[x for x in range(n_s) if split_idx[x]==["2"]]
idxf=[]
idxm=[]
for i in train:
    if parents[i,1] not in idxf:
        idxf=idxf+[parents[i,1]]
    if parents[i,0] not in idxf:
        idxm=idxm+[parents[i,0]]
vidxf=random.sample(idxf,5)
vidxm=random.sample(idxm,5)
#print idxf, idxm, vidxf, vidxm
for i in range(n_s):
    if split_idx[i]==["2"] and ((parents[i,1] in vidxf) or (parents[i,0] in vidxm)):
        print "1"
    else:
        print split_idx[i][0]
