
from scipy.special import wofz # faddeeva func

import numpy as np
class VarPlasmaFunc:

    ##in this class there is Te depenedency. It is not nice because i dont want to inherit plasma class here which is computionally expensive? Ot inherit and iverwrite it here? for now Te = const defined here and dT_dr  = np.array([0,0])

    def __init__(self,Yabs,gamma,local_Nz,local_Nx):

        self.Nz = local_Nz
        self.Nx = local_Nx
        self.Yabs=Yabs
        self.gamma = gamma
    

    

    def zeta(self,n):
        return (1+n*self.Yabs)/self.Nz/self.gamma #should be |Nz| or Nz?

        
    def Z(self,n):
        return 1j*np.pi**0.5*wofz( self.zeta(n) )


    def Zp(self,n):
        #print   1j*np.pi**0.5*( - 2*zeta(r,N,T,n)*wofz(zeta(r,N,T,n)) )# - 2/np.pi**0.5*1j   )
        ## checked!!!
        #print -2*zeta(r,N,T,n)*wofz(zeta(r,N,T,n))
        return 1j*np.pi**0.5*(  -2*self.zeta(n)*wofz(self.zeta(n)) + 2/np.pi**0.5*1j   )

    def Zpp(self,n):
        #print   1j*np.pi**0.5*( - 2*zeta(r,N,T,n)*wofz(zeta(r,N,T,n)) )# - 2/np.pi**0.5*1j   )
        ## checked!!!
        #print -2*zeta(r,N,T,n)*wofz(zeta(r,N,T,n))
        return -2*( self.Z(n)+self.zeta(n)*self.Zp(n) )

    def lam(self):
        return 0.5*self.Nx**2*self.gamma**2/self.Yabs**2
