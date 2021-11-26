from __future__ import division 
from scipy.constants import e,m_e,c
import numpy as np


def X(r):
    z,x=r
    return -2*(x)**2+2


"""
def dX_dr(r):
            
    z,x = r
    dr = 1e-6
    r_dx = np.array([z,x+dr])
    r_dz = np.array([z+dr,x])
    r_dy = np.array([z,x])

    dX_dz=X(r_dz) - X(r)
    dX_dx=X(r_dx) - X(r)
    dX_dy=X(r_dy) - X(r)
        
    return np.array([dX_dz,dX_dx])/dr
"""
def dX_dr(r):
            
    z,x = r
        
    return np.array([np.array([0]),np.array([-4*(x)])])

def gamma(T):
        return np.sqrt(2*T*e/m_e)/c

def dgamma_dr(r):
    z,x = r
    return np.array([0])
