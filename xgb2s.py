#usage: Put the following input files at the same folder as this file:
# geno: space separated snp-index file, first column is number of snp in each interval.
# pheno: 3 columns of phenotype, tab separated
# parents.txt: 2 columns of parental index
# Then do python brr.py geno pheno parents.txt sample_idx 2
# Printout will be the predicted phenotype

import csv
import scipy as SP
import scipy.linalg as LA
import scipy.stats as ST
import xgboost as xgb
import sys
## Main program

n_s = 6210
y = SP.array(list(csv.reader(open(sys.argv[2],'rb'),
                             delimiter='\t'))).astype(float)
y = y[:,int(sys.argv[5])].reshape((n_s,1))#feature 3


# load genotype
X = SP.array(list(csv.reader(open(sys.argv[1],'rb'),delimiter=' '))).astype(float)

n_f = X.shape[0]
X1=X[:,1:]
for i in xrange(n_f):
    X1[i]=(X1[i]-(X1[1]).mean())/X[i,0]
X = X1.T
print X
print X.shape

parents = SP.array(list(csv.reader(open(sys.argv[3],'rb'),
                                   delimiter='\t'))).astype(int)

split_idx=list(csv.reader(open(sys.argv[4],'rb'),delimiter='\t'))
train=[x for x in range(6210) if split_idx[x]==["2"]]
        

idxf=[]
for i in train:
    idxf=idxf+[parents[i,1]]
test1=[]
test2=[]
for i in xrange(n_s):
    if parents[i,1] not in idxf:
        test2=test2+[i]
    elif i not in train:
        test1=test1+[i]
train2=train+test1
test=test1+test2

yhat=SP.zeros((n_s,1))
yhat[train]=y[train]


def train_and_eval(Xtrain,Xtest,ytrain):
    ns=Xtrain.shape[0]
    idx=range(ns)
    SP.random.shuffle(idx)
    train_idx=idx[:int(ns*0.7)]
    valid_idx=idx[int(ns*0.7):]
    xg_train=xgb.DMatrix(Xtrain[train_idx,:],label=ytrain[train_idx,:])
    xg_valid=xgb.DMatrix(Xtrain[valid_idx,:],label=ytrain[valid_idx,:])
    param = {'eta':0.05,'silent':1,
            'subsample':0.5,
           'lambda':0.8}
    model=xgb.train(param, xg_train, 500, [(xg_valid, 'valid')],
                    early_stopping_rounds=30, verbose_eval=10)
    res=model.predict(xgb.DMatrix(Xtest), ntree_limit=model.best_ntree_limit)
    return res.reshape((res.shape[0],1))

yhat[test1]=train_and_eval(X[train], X[test1], yhat[train])
yhat[test2]=train_and_eval(X[train2], X[test2], yhat[train2])


for i in range(len(test)):
    print test[i], yhat[test[i]]
