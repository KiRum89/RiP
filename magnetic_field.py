from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv,jv

from scipy.special import ellipe
from scipy.special import ellipk

class Magnetic_field:
    
    def __init__(self,a,coils):

        self.a = a#0.5 #coils radius
        self.coils = coils#np.array([-0.5,0.5],np.newaxix)



    #what about diff current in coils
    def Yz(self,r):
        z,x = r

        R = np.abs(x/self.a)
        Z = (1/self.a*(self.coils[:,np.newaxis]-z))
        k = ( 4*R/((1+R)**2 + Z**2) )**0.5
        T1 = 1/np.sqrt((1+R)**2 + Z**2)
        T2 = (1-R**2-Z**2)/( (1-R)**2+Z**2 )*ellipe(k**2)+ellipk(k**2)
        return np.sum(T1*T2,axis=0)#*0.45

    def Yx(self,r):
        z,x = r

        R = np.abs(x/self.a)
        
        Z = 1/self.a*(self.coils[:,np.newaxis]-z)
        k = ( 4*R/((1+R)**2 + Z**2) )**0.5
        T1 = Z/R/np.sqrt((1+R)**2 + Z**2)
        #print T1,"R", R,"Z",Z,"x",x
        #print np.sqrt((1+R)**2 + Z**2)
        T2 = (1+R**2+Z**2)/( (1-R)**2+Z**2 )*ellipe(k**2)-ellipk(k**2) # mistake in batchelor?

        
        return -np.sum(T1*T2,axis=0)*np.sign(x)
        """
        if z<0:
        return -np.sum(T1*T2,axis=0)#*0.45
        else:
        return -np.sum(T1*T2,axis=0)#*0.45
        """

    def Yvec(self,r):
        return np.array([self.Yz(r)[0],self.Yx(r)[0]])

    def Yabs(self,r):
        return np.sqrt(self.Yz(r)**2 + self.Yx(r)**2)


    def flux(self,r):
        z,x = r

        R = np.abs(x/self.a)

        Z = 1/self.a*(self.coils[:,np.newaxis]-z)

        k = ( 4*R/((1+R)**2 + Z**2) )**0.5
        T1 = R/np.sqrt((1+R)**2 + Z**2)
        T2 = ((2-k**2)*ellipk(k**2)-2*ellipe(k**2))/k**2
        return np.sum(T1*T2,axis=0)/1.5


    ## do gradient
    def gradY(self,r):
        z,x = r
        dr = 1e-3
        r_dx = np.array([z,x+dr])
        r_dz = np.array([z+dr,x])
        dY_dz=self.Yvec(r_dz) - self.Yvec(r)
        dY_dx=self.Yvec(r_dx) - self.Yvec(r)

        return np.array([dY_dz,dY_dx])/dr


    def dY_dr(self,r):
        return np.dot(self.gradY(r),self.Yvec(r))/self.Yabs(r)
