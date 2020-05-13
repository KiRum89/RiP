from __future__ import division 
from scipy.constants import e,m_e,c
import matplotlib.pyplot as plt
import numpy as np


def X(r):
    p1=2
    p2=3
    X1=0.001
    X0=1.1
    ra = 0.25
    z,x = r
    xb=0.15
    #result = X1 + (X0-X1)*( 1- ( np.sqrt(z**2/5+x**2)/ra )**p1 )**p2
    result = (X0 + (X1-X0)*( 1- (( np.abs(x)-ra )/(xb))**p1 )**p2)
    result[np.where(result>X0)]=X0
    result[np.where( np.abs(x)>ra )]=X1

     
    return result*np.abs(-0.98*z**2 + 1)

def gamma(T):
        return np.sqrt(2*T*e/m_e)/c

def dgamma_dr(r):
    z,x = r
    return np.array([0])

def dX_dr(r):
            
    z,x  = r
    dr = 1e-8
    r_dx = np.array([z,x+dr ])
    r_dz = np.array([z+dr,x ])
    r_dy = np.array([z,x])

    dX_dz=X(r_dz) - X(r)
    dX_dx=X(r_dx) - X(r)

        
    return np.array([dX_dz,dX_dx])/dr

if __name__ == '__main__':
	x = np.linspace(-0.15,0.15,10000)

	plt.plot(x,X([0,x]))
