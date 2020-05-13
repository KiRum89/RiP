import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv, ivp
from scipy.optimize import fsolve

def S(n,lam):
    #print lam#np.exp(-lam)*iv(n,lam)

    return 4*n**2*np.exp(-lam/2)*iv(n,lam/2)/lam # no 2 at the beggining 

def D(q,lam):
    


    n = np.arange(1,10)

    #g = np.exp(-lam)*iv(n,lam)
    #print np.sum( S(n,lam) )
    return 1  -  5*np.sum( S(n,lam)/(q**2-n**2) )

#Y0 = 0.51
#Y1 = 0.99

#Y0 = 0.5001#0.51
#Y1 = 0.50001#0.999



nstep=5000

arr = np.zeros(nstep)


i=0


for l in np.linspace(0.0001,25/2**2,nstep):


    arr[i]=fsolve(D,2.2,l)

    if np.abs(D(arr[i],l))>1e-8:
        print "!!!!warning!!!!"
    i+=1




lam = np.linspace(0.0001,25/2**2,nstep)

plt.plot(np.sqrt(2)*np.sqrt(lam),arr,"o")



def dS_dlam(n,lam):
    return 2*n**2*np.exp(-lam)*(ivp(n,lam)-2*iv(n,lam))/lam**2

def dlam_dw(q,lam):
    return 0#lam/q**2*(0.5*q-1) #remove 1/w

def dq_dw(q):
    return q # = dq_dw * w

def dD_dw(q,lam):
    n = np.arange(1,10)
    return  np.sum( (dS_dlam(n,lam)*dlam_dw(q,lam) - 2*q*dq_dw(q)*S(n,lam))/(q**2-n**2)**2 )  
    
def dD_dk(q,lam):
    n = np.arange(1,10)
    return np.sum( dS_dlam(n,lam)*dlam_dk(q,lam)/(q**2-n**2) ) # remove 1/k
    
def dlam_dk(q,lam):
    return 2*lam # = dlam_dk*k
    
def vg(q,lam):
    return dD_dk(q,lam)/dD_dw(q,lam)


"""
def lam(Y,gamma,Nx):
    return 0.5*Nx**2*gamma**2/Y**2
"""







