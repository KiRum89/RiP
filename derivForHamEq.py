from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from derivOfDispRel import dD_du
from dispRel import disp_rel







def dD_dN(X,Y,gamma,Nz,Nx,dNz_dN,dNx_dN):

    #print "dD_du",dD_du("Nz",X,Y,gamma,Nz,Nx)
    return dD_du("Nx",X,Y,gamma,Nz,Nx)*dNx_dN+dD_du("Nz",X,Y,gamma,Nz,Nx)*dNz_dN



def dD_dr(X,Y,gamma,Nz,Nx,dNz_dr,dNx_dr,dX_dr,dY_dr,dgamma_dr):

    
    
    
    return dD_du("X",X,Y,gamma,Nz,Nx)*dX_dr + dD_du("Y",X,Y,gamma,Nz,Nx)*dY_dr + dD_du("Nx",X,Y,gamma,Nz,Nx)*dNx_dr + dD_du("Nz",X,Y,gamma,Nz,Nx)*dNz_dr + dD_du("gamma",X,Y,gamma,Nz,Nx) * dgamma_dr






#removed 1/omega from the equation because it will cancel in the ray tracing equation
def dD_dw(X,Y,gamma,Nz,Nx):


    return -( 2*X*dD_du("X",X,Y,gamma,Nz,Nx)+Y*dD_du("Y",X,Y,gamma,Nz,Nx)+Nx*dD_du("Nx",X,Y,gamma,Nz,Nx)+Nz*dD_du("Nz",X,Y,gamma,Nz,Nx) )

