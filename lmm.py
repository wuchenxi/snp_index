# gzip -dc mat_snp.gz | python lmm.py

import scipy as SP
import csv
import scipy.linalg as LA
import scipy.optimize as OPT
import sys
from sys import stdin
import numpy as NP

def predict(y_t,K_tt,K_vt,ldelta):
    if y_t.ndim==1:
        y_t = SP.reshape(y_t,(n_train,1))

    n_train = y_t.shape[0]
    delta = SP.exp(ldelta)
    return SP.dot(K_vt,LA.solve(K_tt + delta*SP.eye(n_train),y_t))

def nLLeval(ldelta,Uy,S):
    n_s = Uy.shape[0]
    delta = SP.exp(ldelta)
    
    # evaluate log determinant
    Sd = S+delta
    ldet = SP.sum(SP.log(Sd))

    # evaluate the variance    
    Sdi = 1.0/Sd
    Uy = Uy.flatten()
    ss = 1./n_s * (Uy*Uy*Sdi).sum()

    # evalue the negative log likelihood
    nLL=0.5*(n_s*SP.log(2.0*SP.pi)+ldet+n_s+n_s*SP.log(ss));

    return nLL


def train_lmm(y,K,numintervals=20,ldeltamin=-5,ldeltamax=5,debug=False):
    S,U = LA.eigh(K)
    Uy = SP.dot(U.T,y)
    nllgrid=SP.ones(numintervals+1)*SP.inf
    ldeltagrid=SP.arange(numintervals+1)/(numintervals*1.0)*(ldeltamax-ldeltamin)+ldeltamin
    nllmin=SP.inf
    for i in SP.arange(numintervals+1):
        nllgrid[i]=nLLeval(ldeltagrid[i],Uy,S);
        
    # find minimum
    nll_min = nllgrid.min()
    ldeltaopt_glob = ldeltagrid[nllgrid.argmin()]

    # more accurate search around the minimum of the grid search
    for i in SP.arange(numintervals-1)+1:
        if (nllgrid[i]<nllgrid[i-1] and nllgrid[i]<nllgrid[i+1]):
            ldeltaopt,nllopt,iter,funcalls = OPT.brent(nLLeval,(Uy,S),(ldeltagrid[i-1],ldeltagrid[i],ldeltagrid[i+1]),full_output=True);
            if nllopt<nllmin:
                nllmin=nllopt;
                ldeltaopt_glob=ldeltaopt;
    return ldeltaopt_glob
 
n_s = 6210
y = SP.array(list(csv.reader(open('pheno.csv','r'),
                             delimiter=' '))).astype(float)
y = y[:6210,:]
mid = y[:,0].astype(int)
fid = y[:,1].astype(int)
y = y[:,2].reshape((n_s,1))#feature 1

idxname="./idx/idx"

KL=[]

for row in stdin:
    lrow=[float(x) for x in row.split(' ')]
    KL+=[lrow]
K=SP.array(KL)
K=(4549818.0/SP.trace(K))*K
meig=min(LA.eigh(K)[0])
if meig<0:
    K-=meig*SP.eye(n_s)


for L in range(1,101):

    idxfile=open(idxname+str(L),'r')
    idx1=[int(i) for i in next(idxfile).split(' ')]
    idx2=[int(i) for i in next(idxfile).split(' ')]
    test=[]
    train=[]
    for i in xrange(n_s):
        if SP.isnan(y[i,0]):
            continue
        if fid[i] in idx2 or mid[i] in idx1:
            test=test+[i]
        else:
            train+=[i]
    Kt=K[train][:,train]
    Kv=K[test][:,train]
    yt=4549818.0/(y[train].std())*y[train]
    #ldel=train_lmm(yt, Kt)
    yhat=predict(yt, Kt, Kv, 0)
    corr = 1./len(test) * SP.dot((yhat-yhat.mean()).T,y[test]
                             -y[test].mean())/(yhat.std()*y[test].std())
    print L, corr[0,0]
