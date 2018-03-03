#Usage: python3 summary.py > output_file_name
#For example: python3 summary.py > out
#Input files: emaize_5M.add.gz merbin_infw.txt snpfunc_5M_emaize.gz

import gzip
snp_file_name="emaize_5M.add.gz"
intervals_file_name="merbin_infw.txt"
comments_file_name="snpfunc_5M_emaize.gz"
weights=[1,1,1,1,1,2,3] #weights[0] is the weight for snp that are not classified.

def split_name(s):
    l=str.split(s,".")
    return (int((l[0])[3:]), int((l[1])[2:]))

def int_end(s):
    l=str.split(s.decode("utf-8"),"\t")
    return (int(l[1]),float(l[3]))

snp_f=gzip.open(snp_file_name)
next(snp_f)
interval_f=open(intervals_file_name,"rb")
next(interval_f)
curint=int_end(next(interval_f))
func_f=gzip.open(comments_file_name)
next(func_f)
curfunc=next(func_f)
cursum=[0]*6210

for snp in snp_f:
    snp_l=str.split(snp.decode("utf-8"),",")
    snp_name=split_name(snp_l[0])
    snp_val=[int(x) for x in snp_l[3:]]
    #print(snp_name)
    snp_weight=0
    while True:
        func_l=str.split(curfunc.decode("utf-8"),",")
        if split_name(func_l[0])>=snp_name:
            if split_name(func_l[0])==snp_name:
                #print(func_l[2])
                snp_weight=weights[int(func_l[2])]
            else:
                snp_weight=weights[0]
            break
        curfunc=next(func_f)
    #print(snp_weight)
    if snp_name<curint:
        cursum=[x+snp_weight*y for (x, y) in zip(cursum,snp_val)]
    else:
        for c in cursum[:-1]:
            print(c,end=' ')
        print(cursum[-1])
        while True:
            curint=int_end(next(interval_f))
            if snp_name<curint:
                break
        cursum=[snp_weight*y for y in snp_val]
