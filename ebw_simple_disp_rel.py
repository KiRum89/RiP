import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv, ivp
from scipy.optimize import fsolve

def S(n,lam):
    #print lam#np.exp(-lam)*iv(n,lam)

    return 2*n**2*np.exp(-lam)*iv(n,lam)/lam

def D(Nx,*param):
    
    X,q,gamma = param
    lam = 0.5*Nx**2*gamma**2*q**2
    n = np.arange(1,10)

    #g = np.exp(-lam)*iv(n,lam)

    return -1 +  X*np.sum( S(n,lam)/( 1-(n/q)**2 ) )

#Y0 = 0.51
#Y1 = 0.99

#Y0 = 0.5001#0.51
#Y1 = 0.50001#0.999

q0 = 1.99#1.0101
q1 = 1.001#1.9

nstep=500

arr = np.zeros(nstep)
lam = np.zeros(nstep)

i=0

gamma=0.06256
for q in np.linspace(q0,q1,nstep):
    X =  5*(1/q)**2
    arr[i]=fsolve(D,1.2,args=(X,q,gamma))
    lam[i] = np.sqrt(arr[i]**2*gamma**2*q**2)
    #print D(arr[i],X,Y,1e-2)
    i+=1


plt.plot(lam,np.linspace(q0,q1,nstep),"o")

#lam = np.linspace(0.01,1,1000)

#plt.plot(lam,S(10,lam))



def lam(Y,gamma,Nx):
    return 0.5*Nx**2*gamma**2/Y**2






#plt.plot(arr,np.linspace(0.01,end,nstep))
#plt.plot(np.linspace(0.01,end,nstep),0*np.linspace(0.01,end,nstep))





