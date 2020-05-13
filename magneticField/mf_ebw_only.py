import numpy as np




def Y(z):

    return 0.84033613445+0.001*z**2#0.454318068 + z**2

def q(z):

    return 1/Y(z)

def dq_dz(z):

    return  -2*0.001*z/Y(z)**2




"""
def Yz(self,r):
    z,x = r
    
    R = np.abs(x/self.a)
    Z = (1/self.a*(self.coils[:,np.newaxis]-z))
    k = 
( 4*R/((1+R)**2 + Z**2) )**0.5
    T1 = 1/np.sqrt((1+R)**2 + Z**2)
    T2 = (1-R**2-Z**2)/( (1-R)**2+Z**2 )*ellipe(k**2)+ellipk(k**2)
    return np.sum(T1*T2,axis=0)*0.6

def Yx(self,r):
    z,x = r

    R = np.abs(x/self.a)
    
    Z = 1/self.a*(self.coils[:,np.newaxis]-z)
    k = ( 4*R/((1+R)**2 + Z**2) )**0.5
    T1 = Z/R/np.sqrt((1+R)**2 + Z**2)
    #print T1,"R", R,"Z",Z,"x",x
    #print np.sqrt((1+R)**2 + Z**2)
    T2 = (1-R**2+Z**2)/( (1-R)**2+Z**2 )*ellipe(k**2)-ellipk(k**2)
        
        

         
    if x<0:
        return np.sum(T1*T2,axis=0)*0.6
    else:
        return -np.sum(T1*T2,axis=0)*0.6
            

def Yvec(self,r):
    return np.array([self.Yz(r)[0],self.Yx(r)[0]])

def Yabs(self,r):
    return np.sqrt(self.Yz(r)**2 + self.Yx(r)**2)
    

"""
