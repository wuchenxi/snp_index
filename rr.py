import csv
import scipy as SP
import scipy.linalg as LA
import scipy.stats as ST
import scipy.optimize as OPT

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

def train_nullmodel(y,K,numintervals=100,ldeltamin=-5,ldeltamax=5):
        
    n_s = y.shape[0]

    # rotate data
    S,U = LA.eigh(K)
    Uy = SP.dot(U.T,y)

    # grid search
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

    
    return S,U,ldeltaopt_glob
def predict(y_t,K_tt,K_vt,ldelta):

    n_train = y_t.shape[0]
    
    delta = SP.exp(ldelta)

    print K_vt, y_t

    y_v = SP.dot(K_vt, LA.solve(K_tt + delta*SP.eye(n_train),y_t))
    return y_v


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

K1 = 1.0/n_f*SP.dot(X, X.T)
print K1.shape
K=K1
print K

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

S, U, d=train_nullmodel(yhat[train],K[train][:,train])
delta0=SP.exp(d)
print d
yhat[test1] = predict(y[train],K[train][:,train],K[test1][:,train],d)


S, U, d=train_nullmodel(yhat[train2],K[train2][:,train2])
delta0=SP.exp(d)
print d
yhat[test2] = predict(y[train2],K[train2][:,train2],K[test2][:,train2],d)



print yhat[test]
corr = 1./len(test) * SP.dot((yhat[test]-yhat[test].mean()).T,y[test]
                             -y[test].mean())/(yhat[test].std()*y[test].std())
print corr[0,0]
