import numpy as np
#from various_plasma_functions import VarPlasmaFunc
from scipy.special import iv, ivp # faddeeva func
from scipy.special import wofz # faddeeva func
from plasma_disp_func import g



class HotDiTen:

    def __init__(self,X,Yabs,gamma,local_Nz,local_Nx,num_term):

        self.X = X
        self.Yabs = Yabs
        self.gamma = gamma
        self.Nz = local_Nz
        self.Nx = local_Nx
        self.num_term = num_term
        self.zeta0 = self.zeta(0)
        
        self.n= np.arange(1,num_term+1)
        self.l = self.lam()
        #print self.Yabs, self.Nz
    def Kxx(self):
       
        
        return 1+ self.Bxx()*np.sum(self.n**2*iv(self.n,self.l)*( self.Z(self.n) + self.Z(-self.n) ))

    def Cxy(self):
        return -self.Bxy()*np.sum(self.n*( iv(self.n,self.l) - ivp(self.n,self.l) )*( self.Z(self.n) - self.Z(-self.n) ))
    

    def Cxz(self):
        
        return -self.Bxz()*np.sum(self.n*iv(self.n,self.l)*( -self.Zp(self.n) + self.Zp(-self.n) )/2) #-Zp(q,n) why minus? => because of n
        
    def Kyy(self):
        
        sum=np.sum(( self.n**2*iv(self.n,self.l) + 2*self.l**2*( iv(self.n,self.l)-ivp(self.n,self.l) ) )*(self.Z(self.n)+self.Z(-self.n)))
        return 1 + self.Byy()*( sum + self.l**2*( iv(0,self.l)-ivp(0,self.l) )*(self.Z(0)+self.Z(0)) )

    def Cyz(self):
        
       
        sum=np.sum(( iv(self.n,self.l) - ivp(self.n,self.l) )*( self.zeta(self.n)*self.Z(self.n) + self.zeta(-self.n)*self.Z(-self.n) ))
    
        firstTerm = ( iv(0,self.l) - ivp(0,self.l) )*( self.zeta(0)*self.Z(0) + self.zeta(0)*self.Z(0) )
        return -self.Byz()*(sum + 0.5*firstTerm )

    def Kzz(self):
        
        
        sum=np.sum(iv(self.n,self.l)*( self.zeta(self.n)*self.Zp(self.n) + self.zeta(-self.n)*self.Zp(-self.n) ))
        nn = 0
        firstTerm = iv(nn,self.l)*( self.zeta(nn)*self.Zp(nn) + self.zeta(-nn)*self.Zp(-nn) )
        return 1 - self.Bzz()*( sum + 0.5*firstTerm )
        


    def Bxx(self):

        return  self.X*self.zeta0*np.exp( -self.lam() )/self.lam()

    
    def Bxy(self):

        return self.X*self.zeta0*np.exp(-self.lam())  

    
    def Bxz(self):


        return self.X*np.exp(-self.l)/self.Yabs/self.Nz/self.l  
    

    def Byy(self):

        
        return  self.X*self.zeta0*np.exp(-self.l)/self.l 


    def Byz(self):

        return self.X*np.exp(-self.l)/self.Yabs/self.Nz



    def Bzz(self):


        return self.X*self.zeta0*np.exp(-self.l)  


    



    def zeta(self,n):
        val=(1+n*self.Yabs)/(self.Nz)/self.gamma #should be |Nz| or Nz?
        if val.any() > 1e8:
            print "zeta is very big. Fadeeva wont like it"
        return val

        
    def Z(self,n):
        
        """
        if self.zeta(n)>1e3:
            return -1/self.zeta(n)*(1+1/2/self.zeta(n)**2 ++1/4/self.zeta(n)**4 ) + 1j*np.exp(-self.zeta(n)**2)*np.pi*0.5
        """
        
        #sgn = (-1)**(self.zeta(n)<0)
        sgn = 1
        if self.Nz < 0:
            sgn = -1
        return sgn*1j*np.pi**0.5*wofz(sgn*(self.zeta(n)) )
        

    def Zp(self,n):
        #print   1j*np.pi**0.5*( - 2*zeta(r,N,T,n)*wofz(zeta(r,N,T,n)) )# - 2/np.pi**0.5*1j   )
        ## checked!!!
        #print -2*zeta(r,N,T,n)*wofz(zeta(r,N,T,n))
        #return 1j*np.pi**0.5*(  -2*self.zeta(n)*wofz(self.zeta(n)) + 2/np.pi**0.5*1j   )

    
        sgn = 1
        if self.Nz < 0:
            sgn = -1
      
        return -2 - 2*self.zeta(n)*self.Z(n)#-sgn*1j*np.pi**0.5*2*self.zeta(n)*wofz(np.abs(self.zeta(n))) - 2  

    
    def Zpp_naive(self,n):
        #print   1j*np.pi**0.5*( - 2*zeta(r,N,T,n)*wofz(zeta(r,N,T,n)) )# - 2/np.pi**0.5*1j   )
        ## checked!!!

        return -2*( self.Z(n) + self.zeta(n)*self.Zp(n) )
    
    
    
    def Zpp(self,n):
        # only return the imaginary component
        sgn = 1
        if self.Nz < 0:
            sgn = -1
        
        return sgn*4 * g(sgn*(self.zeta(n))) + 1j*np.imag(self.Zpp_naive(n)) 
    
   

    
    def lam(self):
        
       
        #if val < 1e-7:
        #    print "lambda is very small. Bessel wont like it!"
        #print 0.5*self.Nx**2*self.gamma**2/self.Yabs**2

        return 0.5*self.Nx**2*self.gamma**2/self.Yabs**2









