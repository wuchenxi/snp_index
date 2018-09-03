#Usage: cat Sum_Chrs | python corr.py pheno_6210 sample_idx_va1 2 > Sum
#Here 2 means the third phenotype
#pheno.csv is tab separated, with unknowns as "nan" or "NA"

import fileinput
import sys
import numpy
idx=open(sys.argv[2],"r").readlines()
index=[]
for l in idx:
    index+=[int(l)]
po=open(sys.argv[1],"r").readlines()
npheno=int(sys.argv[3])
pheno=[]
for l in po:
    elm=l.split("\t")[npheno]
    if elm[0]=='N':
        pheno+=[float("nan")]
    else:
        pheno+=[float(elm)]
y=[pheno[i] for i in range(6210) if index[i]>0]
    
count=0
sum=[0.0]*6011
for line in sys.stdin:
    count+=1
    if count==19:
        count=0
        sum=[str(x) for x in sum]
        print(" ".join(sum))
        sum=[0.0]*6211
        continue
    lst=[float(x) for x in line.split(" ")]
    cur_x=[lst[i+1] for i in range(6210) if index[i]>0]
    coef=numpy.corrcoef(cur_x, y)[0,1]
    if coef==coef:
        coef=abs(coef)
        sum=[x+coef*y for x, y in zip(sum, lst)]
    
    
    
