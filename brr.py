#usage: Put the following input files at the same folder as this file:
# geno: space separated snp-index file
# pheno: 3 columns of phenotype, tab separated
# parents.txt: 2 columns of parental index
# Then do python brr.py
# Printout will be the predicted phenotype for samples 4755-6210 from the phenotype of the first 4754 samples

import csv
import scipy as SP
import scipy.linalg as LA
import scipy.stats as ST
import scipy.optimize as OPT
import sklearn.linear_model as LM

## Main program

n_s = 6210
y = SP.array(list(csv.reader(open('pheno','rb'),
                             delimiter='\t'))).astype(float)
y = y[:,2].reshape((n_s,1))#feature 3


# load genotype
X = SP.array(list(csv.reader(open('geno','rb'),delimiter=' '))).astype(float)

n_f = X.shape[0]
for i in xrange(n_f):
    sd=(X[i]).std()
    if sd == 0:
        X[i]=X[i]-(X[i]).mean()
    else:
        X[i]=(X[i]-(X[i]).mean())/sd
X = X.T
print X
print X.shape

parents = SP.array(list(csv.reader(open('parents.txt','rb'),
                                   delimiter='\t'))).astype(int)

idxf=[]
train=list(range(4754))
for i in train:
    idxf=parents[i,m]
test1=[]
test2=[]
for i in xrange(n_s):
    if parents[i,1] in idxf:
        test2=test2+[i]
    else:
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


for i in range(len(test)):
    print test[i], yhat[test[i]]
#corr = 1./len(test) * SP.dot((yhat[test]-yhat[test].mean()).T,y[test]
#                             -y[test].mean())/(yhat[test].std()*y[test].std())
#print corr[0,0]
