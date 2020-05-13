from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import iv
from plasma import PlasmaParabMF, SlabPlasma


s22 = []
s11 = []
tt = []

n = np.arange(-2,2)

def lam(q,beta,nPerp):
    
    return 0.5*(nPerp*q*beta)**2

def theta(X,q,l):

    return X*np.exp(-l)*iv(n,l)*q**2/l

def nParall2(q,l,beta,theta,X):

    t = 1 + q**2*X/l#*np.sum(iv(n,l))*np.exp(-l)
    s1 = np.sum(q*theta/(q-n))
    s2 = np.sum(theta*q**3*beta**2/2/(q-n)**3)
    s22.append( ( 1/(q-n) ) )

    s11.append(s1)
    tt.append(t)
    return  (t - s1)/s2 

def xi(q,beta,nParall):

    return ((q - n)/q/beta/nParall)

def simpDispRel(q,l,X,beta,nParall):
    return 1 + X*q**2/l - np.sum( theta(X,q,l) * ( 1/xi(q,beta,nParall)  + 1/2/xi(q,beta,nParall)**3 )/beta/nParall ) 

def simpDispRelEBW(q,l,X):
   
    print   X*q**2/l - np.sum(theta(X,q,l) * ( q**2/(q**2-n**2)) )
    return 1 + X*q**2/l - np.sum( theta(X,q,l) * ( q**2/(q**2-n**2)) )

#s = np.load("paraMF.npy")
s = np.load("volpeDispRelPlot.npy")
i = 0
p = SlabPlasma(1)#PlasmaParabMF(1)
nparal_plt = []# np.zeros((1000,1000))
eps = []
ll = []

"""
for i in  range(0,s[:,0].size):
    z,x = s[i,0],s[i,1]
    r = [z,x]
    q = 1/p.Yabs(r)
    b = p.gamma(r)
    X = p.X(r)
    l = lam(q,b,s[i,3]-0.2)
    th = theta(X,q,l)
    xii = xi(q,b,s[i,2])
    nparal_plt.append(nParall2(q,l,b,th,X))
    #nParall = s[i,2] 
    #eps.append( simpDispRel(q,l,xii,X,b,s[i,2]))#np.abs(nParall) ))
    ll.append(l)



"""

for nx in np.linspace(10,99,1000):



    b = 0.01
    X = 2 #1.4
    q= 1/0.68 #1.1
    l = lam(q,b,nx)
    th = theta(X,q,l)

    nparal_plt.append(nParall2(q,l,b,th,X))
    ll.append(nx)





"""
b = 0.04
X = 2
i=0
j=0

for q in np.linspace(1.01,1.9,1000):
    # = 0
    #or nPerp in np.linspace(1,50,1000):

        

    l = lam(q,b, 40)
    th = theta(X,q,l)
    nparal_plt.append(nParall2(q,l,b,th,X)**0.5)
    ll.append(q)
    
    # +=1
        
    # += 1

l = lam(1.19938,0.04,40)

print "l",l,simpDispRel(1.19938,l,2,0.04,0.059098754815139855)#0.028)#0.08877651509632188)

"""
