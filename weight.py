#Usage: cat Sum_Chrs | python corr.py > Sum
#Here 2 means the third phenotype
#pheno.csv is tab separated, with unknowns as "nan" or "NA"

import fileinput
import sys
import numpy
weights=[0,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]
#0 followed by the chosen weights of 18 classes    
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
    sum=[x+weights[count]*y for x, y in zip(sum, lst)]
    
    
    
