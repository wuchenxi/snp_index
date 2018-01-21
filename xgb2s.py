import csv
import scipy as SP
import scipy.linalg as LA
import scipy.stats as ST
import xgboost as xgb

## Main program

n_s = 4754
y = SP.array(list(csv.reader(open('pheno','rb'),
                             delimiter='\t'))).astype(float)
y = y[:4754,:]
y = y[:,2].reshape((n_s,1))#feature 3


# load genotypes
X = SP.array(list(csv.reader(open('geno','rb'),delimiter='\t'))).astype(float)

# remove snp label
X = X[:,:n_s]
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

#use sample_idx
split_idx=SP.array(list(csv.reader(open('sample_idx','rb'),
                                   delimiter='\t'))).astype(int)
idxmc=[]
idxfc=[]
for i in xrange(split_idx.shape[0]):
    if split_idx[i,0]==2:
        print i, parents[i,1]
        if parents[i,1] not in idxfc:
            idxfc=idxfc+[parents[i,1]]
        if parents[i,0] not in idxmc:
            idxmc=idxmc+[parents[i,0]]
idxm=[x for x in range(1,191) if x not in idxmc]
idxf=[x for x in range(1,26) if x not in idxfc]
print idxm, idxf

#random index
##parents=parents[:4754,:]
##idxm=range(1,191)
##SP.random.shuffle(idxm)
##idxm=idxm[:5]
##idxf=range(1, 26)
##SP.random.shuffle(idxf)
##idxf=idxf[:5]

train=[]
test1=[]
test2=[]
for i in xrange(n_s):
    if parents[i,1] in idxf:
        test2=test2+[i]
    else:
        if parents[i,0] in idxm:
            test1=test1+[i]
        else:
            train=train+[i]
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
    model=xgb.train(param, xg_train, 300, [(xg_valid, 'valid')],
                    early_stopping_rounds=30, verbose_eval=10)
    res=model.predict(xgb.DMatrix(Xtest), ntree_limit=model.best_ntree_limit)
    return res.reshape((res.shape[0],1))


yhat[test1]=train_and_eval(X[train], X[test1], yhat[train])
yhat[test2]=train_and_eval(X[train2], X[test2], yhat[train2])


print yhat[test].ravel()
corr = 1./len(test) * SP.dot((yhat[test]-yhat[test].mean()).T,y[test]
                             -y[test].mean())/(yhat[test].std()*y[test].std())
print corr[0,0]
