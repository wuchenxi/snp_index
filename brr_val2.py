# Usage is the same as brr_val.py, output a table of 3 columns:
# column 1: index of samples in the validation set.
# column 2: predicted value in the validation set.
# column 3: true value in the validation set.

import csv
import scipy as SP
import scipy.linalg as LA
import scipy.stats as ST
import scipy.optimize as OPT
import sklearn.linear_model as LM
import sys
## Main program

n_s = 6210
y = SP.array(list(csv.reader(open(sys.argv[2],'r'),
                             delimiter='\t'))).astype(float)
y = y[:,int(sys.argv[5])].reshape((n_s,1))#feature 3


# load genotype
X = SP.array(list(csv.reader(open(sys.argv[1],'r'),delimiter=' '))).astype(float)

n_f = X.shape[0]
X1=X[:,1:]
for i in range(n_f):
    X1[i]=(X1[i]-(X1[1]).mean())/X[i,0]
X = X1.T
#print X
#print X.shape

parents = SP.array(list(csv.reader(open(sys.argv[3],'r'),
                                   delimiter='\t'))).astype(int)

split_idx=list(csv.reader(open(sys.argv[4],'r'),delimiter='\t'))
train=[x for x in range(6210) if split_idx[x]==["2"]]
valid=[x for x in range(6210) if split_idx[x]==["1"]]

idxf=[]
for i in train:
    if parents[i,1] not in idxf:
        idxf=idxf+[parents[i,1]]
test1=[]
test2=[]
for i in range(n_s):
    if parents[i,1] not in idxf:
        test2=test2+[i]
    elif i not in train:
        test1=test1+[i]
train2=train+test1
test=test1+test2

yhat=SP.zeros((n_s,1))
yhat[train]=y[train]

def train_and_eval(Xtrain,Xtest,ytrain):
    ytrain=ytrain.ravel()
    reg=LM.BayesianRidge(n_iter=100)
    reg.fit(Xtrain, ytrain)
    res=reg.predict(Xtest)
    return res.reshape((res.shape[0],1))


yhat[test1]=train_and_eval(X[train], X[test1], yhat[train])
yhat[test2]=train_and_eval(X[train2], X[test2], yhat[train2])

#summ=[(y[test[i]]-yhat[test[i]])**2 for i in range(len(test)) if test[i] in valid]
for s in valid:
    print(s, yhat[s][0], y[s][0])
#print sum(summ)
