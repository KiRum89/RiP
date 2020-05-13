from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from coldDerivOfDispRel import dD_du






def dD_dN(X,Y,Nz,Nx,dNz_dN,dNx_dN,m):
    return dD_du("Nx",X,Y,Nz,Nx,m)*dNx_dN+dD_du("Nz",X,Y,Nz,Nx,m)*dNz_dN




def dD_dr(X,Y,Nz,Nx,dNz_dr,dNx_dr,dX_dr,dY_dr,m):


    
    return dD_du("X",X,Y,Nz,Nx,m)*dX_dr + dD_du("Y",X,Y,Nz,Nx,m)*dY_dr + dD_du("Nx",X,Y,Nz,Nx,m)*dNx_dr + dD_du("Nz",X,Y,Nz,Nx,m)*dNz_dr






#removed 1/omega from the equation because it will cancel in the ray tracing equation
def dD_dw(X,Y,Nz,Nx,m):

    return -( 2*X*dD_du("X",X,Y,Nz,Nx,m)+Y*dD_du("Y",X,Y,Nz,Nx,m)+Nx*dD_du("Nx",X,Y,Nz,Nx,m)+Nz*dD_du("Nz",X,Y,Nz,Nx,m) )

