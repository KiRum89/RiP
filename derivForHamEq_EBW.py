from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
#from derivOfDispRel import dD_du
from dispRel_EBW import disp_rel
#import dispRel_EBW



def dD_dX(X,Y,gamma,Nz,Nx):
    dX = 1e-6
    return (disp_rel(X+dX,Y,gamma,Nz,Nx)-disp_rel(X,Y,gamma,Nz,Nx))/dX

def dD_dY(X,Y,gamma,Nz,Nx):
    dY = 1e-6
    return (disp_rel(X,Y+dY,gamma,Nz,Nx)-disp_rel(X,Y,gamma,Nz,Nx))/dY

def dD_dNz(X,Y,gamma,Nz,Nx):
    dNz = 1e-6
    return (disp_rel(X,Y,gamma,Nz+dNz,Nx)-disp_rel(X,Y,gamma,Nz,Nx))/dNz

def dD_dNx(X,Y,gamma,Nz,Nx):
    dNx = 1e-6
    return (disp_rel(X,Y,gamma,Nz,Nx+dNx)-disp_rel(X,Y,gamma,Nz,Nx))/dNx

def dD_dgamma(X,Y,gamma,Nz,Nx):
    dgamma = 1e-6
    return (disp_rel(X,Y,gamma+dgamma,Nz,Nx)-disp_rel(X,Y,gamma,Nz,Nx))/dgamma




def dD_dN(X,Y,gamma,Nz,Nx,dNz_dN,dNx_dN):
    return dD_dNx(X,Y,gamma,Nz,Nx)*dNx_dN+dD_dNz(X,Y,gamma,Nz,Nx)*dNz_dN



def dD_dr(X,Y,gamma,Nz,Nx,dNz_dr,dNx_dr,dX_dr,dY_dr,dgamma_dr):

    
    
    
    return dD_dX(X,Y,gamma,Nz,Nx)*dX_dr + dD_dY(X,Y,gamma,Nz,Nx)*dY_dr + dD_dNx(X,Y,gamma,Nz,Nx)*dNx_dr + dD_dNz(X,Y,gamma,Nz,Nx)*dNz_dr + dD_dgamma(X,Y,gamma,Nz,Nx) * dgamma_dr






#removed 1/omega from the equation because it will cancel in the ray tracing equation
def dD_dw(X,Y,gamma,Nz,Nx):


    return -( 2*X*dD_dX(X,Y,gamma,Nz,Nx)+Y*dD_dY(X,Y,gamma,Nz,Nx)+Nx*dD_dNx(X,Y,gamma,Nz,Nx)+Nz*dD_dNz(X,Y,gamma,Nz,Nx) )
