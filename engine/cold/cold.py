import sys
import os
sys.path.insert(0,os.path.abspath(""))

import refractive_index as ri
from derivForHamEq_c import *

import conf
from conf import mf,p
from  scipy.integrate import odeint
from coldDispRel import cold_disp_rel,coldNx_Xm,coldDispRel_pol,S,coldNx_Om
import params



Num=1
m = params.m
def get_ray(init_cond,ncycle):

    #init_cond=np.array([z0,x0,nz0,nx0])   
    def f(q,t):
        z,x,Nz,Nx = q.reshape(4,Num)
        # all the variables are in global CS
        N=np.array([Nz,Nx])
        r = np.array([z,x])
        X=p.X(r)
        Y=mf.Yabs(r)

        dnz_dN=ri.dnz_dN(r,N)
        dnx_dN=ri.dnx_dN(r,N)
        dnz_dr=ri.dnz_dr(r,N)
        dnx_dr=ri.dnx_dr(r,N)
        dX_dr = p.dX_dr(r)
        dY_dr = mf.dY_dr(r)
       
        nz,nx=ri.nz(r,N),ri.nx(r,N)
      
        n = np.array([nz,nx])

        
        vg = -np.real( dD_dN( X,Y,nz,nx,dnz_dN,dnx_dN,m  ) / dD_dw( X,Y,nz,nx,m  ) )
        
        a = np.real( dD_dr(X,Y,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,m) / dD_dw( X,Y,nz,nx,m  ) )
        arr = np.array([vg,a])

        print(cold_disp_rel(X,Y,nz,nx,m), vg[0],vg[1])

        return arr.reshape(4*Num)


    def solv(init_cond,t):
        # initial cond
    
      
        init = init_cond
        print(init)
        soln = odeint(f,np.reshape(init,4*Num),t,full_output=0, printmessg=1) 


        return soln
    #ncycle = 5
    nstep = 50000 #number of steps is important. Too much can be bad apparently due too floating ariphmetics. 2e4 didnt work, 20 work well???

    

    t  = np.linspace(0, ncycle, nstep)
    soln=solv(init_cond,t)
    return soln


if __name__=='__main__':
	Num = 1
	m="Xm"
	z0 = np.array([0.01])
	x0 = np.array([-0.11])#np.array([1.01])
	r = np.array([z0,x0])
	


	Nz0 = np.array([0.01])#np.array([ (1/np.sqrt(2) / (1 + 1/np.sqrt(2)))**0.5 - 0.1 ] )

	ncycle =1 


	Y=mf.Yabs(r)
	X = p.X(r)
	Nx0 = coldNx_Xm(X,Y,Nz0)
	N = np.array([Nz0,Nx0])
	nz,nx=ri.nz(r,N),ri.nx(r,N)
	
	init_cond=np.array([z0,x0,Nz0,Nx0])  
	soln=get_ray(init_cond,ncycle)	
