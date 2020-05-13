from  __future__ import division
import numpy as np

from selector import H


def dD_du(u,X,Yabs,Nz,Nx,m):
    sin_theta=Nx/np.sqrt(Nx**2+Nz**2)

    return -2*( Nx*H(u,"Nx")+Nz*H(u,"Nz")  +  (  (1-2*X)*Q(X,Yabs,Nz,Nx,m)*H(u,"X")-X*(1-X)*dQ_du(u,X,Yabs,Nz,Nx,m)  )/Q(X,Yabs,Nz,Nx,m)**2  )

def sinTh(Nz,Nx):
    #
    return Nx/np.sqrt( Nx**2 + Nz**2 )

def cosTh(Nz,Nx):
    
    return Nz/np.sqrt( Nx**2 + Nz**2 )



def G(X,Yabs,Nz,Nx):
    return np.sqrt( Yabs**4*sinTh(Nz,Nx)**4 + 4*(1-X)**2*Yabs**2*cosTh(Nz,Nx)**2 )



def Q(X,Yabs,Nz,Nx,m):
    
    if m=="Om":
        return 2*(1-X)-Yabs**2*sinTh(Nz,Nx)**2+G(X,Yabs,Nz,Nx)
    if m=="Xm" :
        return 2*(1-X)-Yabs**2*sinTh(Nz,Nx)**2-G(X,Yabs,Nz,Nx)


def disp_rel(X,Yabs,Nz,Nx):
    return 1-Nx**2-Nz**2 - 2*X*(1-X)/Q(X,Yabs,Nz,Nx,m)


def dQ_du(u,X,Yabs,Nz,Nx,m):
    
    
    if m=="Om":
        return   -( 2*H(u,"X")+2*Yabs*sinTh(Nz,Nx)**2*H(u,"Y")+Yabs**2*dSin2Th_du(u,Nz,Nx) - dG_du(u,X,Yabs,Nz,Nx) )
    else: return -( 2*H(u,"X")+2*Yabs*sinTh(Nz,Nx)**2*H(u,"Y")+Yabs**2*dSin2Th_du(u,Nz,Nx) + dG_du(u,X,Yabs,Nz,Nx) )


def dSin2Th_du(u,Nz,Nx):
    return 2*Nx*Nz*( Nz*H(u,"Nx")-Nx*H(u,"Nz") ) / ( Nx**2+Nz**2 )**2

def dG_du(u,X,Yabs,Nz,Nx):
    
    return ( 2*Yabs**3*sinTh(Nz,Nx)**4*H(u,"Y") + Yabs**4*sinTh(Nz,Nx)**2*dSin2Th_du(u,Nz,Nx)\
             + 2*Yabs*(1-X) * ( 2*cosTh(Nz,Nx)**2*H(u,"Y")*(1-X)\
                                - Yabs*(  2*H(u,"X")*cosTh(Nz,Nx)**2\
                                          + (1-X)*dSin2Th_du(u,Nz,Nx)) ) )/G(X,Yabs,Nz,Nx)
