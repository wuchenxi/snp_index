import csv
import scipy as SP
import scipy.linalg as LA
import scipy.stats as ST
import scipy.optimize as OPT
import sklearn.linear_model as LM

## Main program

n_s = 4754
y = SP.array(list(csv.reader(open('pheno','rb'),
                             delimiter='\t'))).astype(float)
y = y[:4754,:]
y = y[:,2].reshape((n_s,1))


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
parents=parents[:4754,:]
idxm=range(1,191)
SP.random.shuffle(idxm)
idxm=idxm[:5]
idxf=range(1, 26)
SP.random.shuffle(idxf)
idxf=idxf[:5]

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

reg=LM.BayesianRidge(n_iter=20)
reg.fit(X[train],yhat[train][:,0])
res=reg.predict(X[test1])
yhat[test1]=res.reshape((res.shape[0],1))

reg=LM.BayesianRidge(n_iter=20)
reg.fit(X[train2],yhat[train2][:,0])
res=reg.predict(X[test2])
yhat[test2]=res.reshape((res.shape[0],1))


print yhat[test]
corr = 1./len(test) * SP.dot((yhat[test]-yhat[test].mean()).T,y[test]
                             -y[test].mean())/(yhat[test].std()*y[test].std())
print corr[0,0]
