#Usage: python summary.py > output_file_name
#For example: python summary.py > out
#Input files: emaize_5M.add.gz merbin_infw.txt snpfunc_5M_emaize.gz 200kpos.txt
#Then "out" is a space separated snp-index file and we can use various ML methods to analyze it.

import gzip
import sys
snp_file_name="emaize_5M.add.gz"
intervals_file_name="merbin_infw.txt"
comments_file_name="snpfunc_5M_emaize.gz"
exclude_file_name="200kpos.txt"
weights=[1,1,1,1,2,3] 

def split_name(s):
    l=s.split(".")
    return (int((l[0])[3:]), int((l[1])[2:]))

def int_end(s):
    l=s.decode("ascii").split("\t")
    return (int(l[1]),float(l[3]))
def ind_exclude(s):
    l=s.decode("ascii").split("\t")
    return (int(l[0]),float(l[1]))

snp_f=gzip.open(snp_file_name)
next(snp_f)
interval_f=open(intervals_file_name,"rb")
next(interval_f)
curint=int_end(next(interval_f))
func_f=gzip.open(comments_file_name)
next(func_f)
func_l=(next(func_f)).decode("ascii").split(",")
cursum=[0]*6210

exclude_f=open(exclude_file_name,"rb")
curex=ind_exclude(next(exclude_f))

for snp in snp_f:
    snp_l=snp.decode("ascii").split(",")
    snp_name=split_name(snp_l[0])
    if snp_name==curex:
        continue
    while snp_name>curex:
        exclude_s=exclude_f.readline()
        if exclude_s==b'':
            curex=(11,1)
        else:
            curex=ind_exclude(exclude_s)
    snp_val=[int(x) for x in snp_l[3:]]
    #print(snp_name)
    snp_weight=0
    while split_name(func_l[0])<snp_name:
        func_s=func_f.readline()
        if func_s==b'':
            func_l=["chr11.s_1",0,1]
        else:
            func_l=(func_s).decode("ascii").split(",")
    if split_name(func_l[0])==snp_name:
        snp_weight=weights[int(func_l[2])-1]
    else:
        continue
    if snp_name<curint:
        cursum=[x+snp_weight*y for (x, y) in zip(cursum,snp_val)]
    else:
        print(' '.join(map(str,cursum)))
        while curint<=snp_name:
            interval_s=next(interval_f)
            if interval_s==b'':
                sys.exit(0)
            else:
                curint=int_end(interval_s)
        cursum=[snp_weight*y for y in snp_val]
