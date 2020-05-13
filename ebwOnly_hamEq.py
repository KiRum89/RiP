from __future__ import division
import numpy as np
from scipy.special import iv,wofz
from conf import mf

def Z(xi):
    return 1j*np.pi**0.5*wofz(xi)

def Zp(xi):
    return -2 - 2*xi*Z(xi)

def lam(q,gamma,nPerp):
    return  0.5*(q*gamma*nPerp)**2

def theta(n,X,q,gamma,nPerp):
    l = lam(q,gamma,nPerp)
    return 2*X*np.exp(-l)*iv(n,l)/(nPerp*gamma)**2

def eps(X,q,gamma,nParal, nPerp):
    
    n = np.arange(-50,50)
    xi = (q-n)/gamma/q/nParal
    th=theta(n,X,q,gamma,nPerp)

    #return 1 + 2*X/(nPerp*gamma)**2 + np.sum(th*(Z(xi))/gamma/nParal)
    return 1 + 2*X/(nPerp*gamma)**2 - np.sum(th*( q/(q-n) + 0.5*q**3/(q-n)**3*(gamma*nParal)**2 ))

##derivatives
def wdeps_dw(X,q,gamma,nParal, nPerp):

    d = 1e-8
    #these derivatives are multiplied by w to get rid off the 1/w multiplier
    dX_dw = - 2*X
    dq_dw =  q

    


    return (eps(X+d,q,gamma,nParal,nPerp) - eps(X,q,gamma,nParal, nPerp) )/d * dX_dw\
        +( eps(X,q+d,gamma,nParal,nPerp) - eps(X,q,gamma,nParal, nPerp) )/d * dq_dw\
        -( eps(X,q,gamma,nParal+d,nPerp) - eps(X,q,gamma,nParal, nPerp) )/d * nParal\
        -( eps(X,q,gamma,nParal,nPerp+d) - eps(X,q,gamma,nParal, nPerp) )/d * nPerp
    

def deps_dz(X,q,gamma,nParal, nPerp,z):
    d = 1e-8
    dq_dz = mf.dq_dz(z)
    #print "mf",dq_dz
    return ( eps(X,q+d,gamma,nParal,nPerp) - eps(X,q,gamma,nParal, nPerp) )/d * dq_dz


def deps_dnParal(X,q,gamma,nParal, nPerp,z):
    d = 1e-8
    #return ( eps(X,q,gamma,nParal+d,nPerp) - eps(X,q,gamma,nParal, nPerp) )/d

    n = np.arange(-50,50)
    #xi = (q-n)/gamma/q/nParal
    #print "xi",xi
    th=theta(n,X,q,gamma,nPerp)
    #dxi_dnParal = -(q - n)/gamma/nParal**2/q
    return ( eps(X,q,gamma,nParal+d,nPerp) - eps(X,q,gamma,nParal, nPerp) )/d

    #return np.sum(th*Zp(xi)*dxi_dnParal-gamma*Z(xi))/(gamma*nParal)**2





