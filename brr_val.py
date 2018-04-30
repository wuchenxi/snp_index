#usage: Put the following input files at the same folder as this file:
# geno: space separated snp-index file, first column is number of snp in each interval.
# pheno: 3 columns of phenotype, tab separated
# parents.txt: 2 columns of parental index
# sample_idx: Index file indicating the train/test/validation split
# Then do python brr_val.py geno pheno parents.txt sample_idx 2
# The last number in the printout will be the covariance of predicted value on validation set

import csv
import scipy as SP
import scipy.linalg as LA
import scipy.stats as ST
import scipy.optimize as OPT
import sklearn.linear_model as LM
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
valid=[x for x in range(6210) if split_idx[x]==["1"]]

idxf=[]
for i in train:
    if parents[i,1] not in idxf:
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

def train_and_eval_with_select(Xtrain,Xtest,ytrain):
    ytrain=ytrain.ravel()
    from sklearn.feature_selection import f_regression
    F, p=f_regression(Xtrain,ytrain)
    features=[]
    for i in xrange(Xtrain.shape[1]):
        if p[i]<0.001:
            print p[i]
            features+=[i]
    print len(features)
    Xtrain=Xtrain[:,features]
    Xtest=Xtest[:,features]
    reg=LM.BayesianRidge(n_iter=30)
    reg.fit(Xtrain, ytrain)
    res=reg.predict(Xtest)
    return res.reshape((res.shape[0],1))

def train_and_eval(Xtrain,Xtest,ytrain):
    ytrain=ytrain.ravel()
    reg=LM.BayesianRidge(n_iter=100)
    reg.fit(Xtrain, ytrain)
    res=reg.predict(Xtest)
    return res.reshape((res.shape[0],1))


yhat[test1]=train_and_eval(X[train], X[test1], yhat[train])
yhat[test2]=train_and_eval(X[train2], X[test2], yhat[train2])

summ=[(y[test[i]]-yhat[test[i]])**2 for i in range(len(test)) if test[i] in valid]
print sum(summ)
