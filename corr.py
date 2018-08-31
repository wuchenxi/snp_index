#Usage: gzip -dc sum_chr1.gz | python3 corr.py pheno.csv index 2 | gzip > chr1.gz
#Here 2 means the third phenotype
#pheno.csv is comma separated, with unknowns as "nan"

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
    pheno+=[float(l.split(",")[npheno])]
y=[pheno[i] for i in range(6210) if index[i]>0]
    
count=0
sum=[0.0]*6011
for line in sys.stdin.readlines():
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
    
    
    
