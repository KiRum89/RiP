import numpy as np


class Local_N:

    def __init__(self,Yvec,gradY, global_N):
        self.Yvec = Yvec
        self.Yabs = np.sqrt(np.sum(Yvec**2))
        self.gradY = gradY
        self.global_N = global_N

        
    def Nz(self):

        if self.Nz < 1e-6:
            print "Nz is small"

        return np.dot(self.global_N,self.Yvec/self.Yabs)

    def Nx(self):

        ## raise exception if Nx is too small

        if self.Nx < 1e-5:
            print self.Nx, "Nx is small"
        return np.sqrt(np.sum(self.global_N**2)-self.Nz()**2)

    
    def dNz_dN(self):


        return self.Yvec/self.Yabs

    def dNx_dN(self):
        return ( self.global_N - self.Yvec/self.Yabs*self.Nz() )/self.Nx()

    def dNz_dr(self):
        #print self.gradY,self.Nx()/self.Yabs*self.dNx_dN()
        return np.dot( self.gradY,self.Nx()/self.Yabs*self.dNx_dN() ) #np.array([0,0])

    def dNx_dr(self):
        
        #print -self.Nz()/self.Nx(),self.dNz_dr()
       
        return -self.Nz()/self.Nx()*self.dNz_dr()
