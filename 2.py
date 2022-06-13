# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 11:40:40 2021

@author: duduv
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 15:27:17 2021

@author: vincent.d
"""

import os
import numpy as np
import csv
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.signal import savgol_filter
from matplotlib import pyplot as plt

from scipy.sparse import csc_matrix, eye, diags





def baseline_als_optimized(y, lam, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just *update diagonal values
        Z = W + D
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

def airPLS(x, lambda_=100, porder=1, itermax=15):
    m=x.shape[0]
    #print(m)
    w=np.ones(m)
    for i in range(1,itermax+1):
        z=WhittakerSmooth(x,w,lambda_, porder)
        d=x-z
        dssn=np.abs(d[d<0].sum())
        if(dssn<0.001*(abs(x)).sum() or i==itermax):
            if(i==itermax): print('WARING max iteration reached!')
            break
        w[d>=0]=0 # d>0 means that this point is part of a peak, so its weight is set to 0 in order to ignore it
        w[d<0]=np.exp(i*np.abs(d[d<0])/dssn)
        w[0]=np.exp(i*(d[d<0]).max()/dssn) 
        w[-1]=w[0]
    return z

def WhittakerSmooth(x,w,lambda_,differences=1):
    X=np.matrix(x)
    m=X.size
    i=np.arange(0,m)
    E=eye(m,format='csc')
    D=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
    W=diags(w,0,shape=(m,m))
    A=csc_matrix(W+(lambda_*D.T*D))
    B=csc_matrix(W*X.T)
    background=spsolve(A,B)
    return np.array(background)





mot_cles='.txt'
arr1 = os.listdir('Spectres_a_tracer')
print(arr1)





Xmin=0
Xmax=3000

for a in arr1:
    if '.csv' in a : ##Boucle si csv
        with open('Spectres_a_tracer/'+str(a)) as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            data1=[row for row in reader]
            x0=data1[21][:]
            y0=data1[23][:]
        del y0[-1]
        x1=np.array(x0,dtype='float')
        y1=np.array(y0,dtype='float')
    if '.txt' in a : ##Boucle si csv
        a1=np.loadtxt('Spectres_a_tracer/'+str(a))
        x=a1[0,:]
        y=a1[1,:]
        x1=[]
        y1=[]
        for i in range(0,np.size(y)-1):
           if x[i]>Xmin and x[i]<Xmax :
            x1=x1+[x[i]]
            y1=y1+[y[i]]
    
    ###Boucle de mofidication si jamai
    #print(x1)
    x1=np.array(x1)
    y1=np.array(y1)
    # ytest=airPLS(y1,lambda_=50, porder=3, itermax=100) 
    # ytest=y1-ytest
    # y_corr=baseline_als_optimized(ytest,532,0.00000009)###Pas mal du tout 
    # y2=ytest-y_corr
     
 
      
    # #yhat =savgol_filter(ytest, 15, 8)
    # a=a.replace('.txt','')
    
    # yhat=y2
    ###☻Normalisation a utiliser
    #yhat=yhat/yhat[333]
    #yhat=yhat/yhat[690]
    
    
    #yhat=yhat/y[2127] ####•Meilleur correction pour sorbitol/dextrose Ycorr VD 28 10 2021 avec Y_corr
    #yhat=yhat/y[2127]
    y1=y1/y1[1065]
    
    y1=y1/y1[2531]
    
    #yhat=yhat/yhat[2745]
    
    
    
    
    #yhat=yhat/np.max(yhat)
    #yhat=yhat/y[2134]
    plt.plot(x1,y1,label=a)
    
    plt.legend()

# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.xlabel('Raman shift(cm$^{-1}$)',fontsize=25)
# plt.ylabel('Normalized Raman Intensity (a.u) ',fontsize=25)
plt.show()
      
 
        
