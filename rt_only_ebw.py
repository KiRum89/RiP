from ebwOnly_hamEq import *
from mf_ebw_only import q
from  scipy.integrate import odeint
#import numpy as np
import matplotlib.pyplot as plt 

z0  = 0
nParal0 =0.059098754815139855# 0.088776515096321881#0.08878
nPerp = 40
gamma = 0.04
X = 2


init_cond = [z0,nParal0]


vgz = []

def check_xi(q,nParal):
        n = np.arange(-50,50)
        xi = (q-n)/gamma/q/nParal
        # if xi[np.where(np.abs(xi)<1)[0][0]]:
        print xi[np.where(xi<1)[0][0]]#"wanring"
        

def fX(r):
        z,x= r
        return np.abs(x)


def get_ray(init_cond):

    def f(a,tau):
    
        z = a[0]
        nParal = a[1]

        check_xi(q(z),nParal)
        #print eps(X,q(z),gamma,nParal,nPerp)
        f0 = -np.real(deps_dnParal(X,q(z),gamma,nParal,nPerp,z)/wdeps_dw(X,q(z),gamma,nParal,nPerp))

        
        f1 = np.real(deps_dz(X,q(z),gamma,nParal,nPerp,z)/wdeps_dw(X,q(z),gamma,nParal,nPerp))
        
        
        vgz.append(f0) 

        return [f0,f1]


    def solv(init_cond,t):
        # initial cond
    
      
        init = init_cond
        print init
        soln = odeint(f,init,t,full_output=0, printmessg=1) 


        return soln
    ncycle = 1000
    nstep = 400000 #number of steps is important. Too much can be bad apparently due too floating ariphmetics. 2e4 didnt work, 20 work well???

    

    t  = np.linspace(0, ncycle, nstep)
    soln=solv(init_cond,t)
    return soln


soln=get_ray(init_cond)




z = 0
nParal = nParal0
nPerp = 40

print  "l",lam(q(z),0.04, nPerp), "eps",eps(X,q(z),gamma,nParal, nPerp)

print "vgz",deps_dnParal(X,q(z),gamma,nParal,nPerp,z),wdeps_dw(X,q(z),gamma,nParal,nPerp)

