import sys
import os
##set paths to packages
project = "RiP_general_new"
sys.path.insert(0,os.path.abspath(""))
import conf
from  scipy.integrate import odeint



import params

import refractive_index as ri
from conf import p,mf


from derivForHamEq import *
from dispRel import A,B,C, disp_rel
from coldDispRel import cold_disp_rel,coldNx_Xm,coldDispRel_pol,S,coldNx_Om
from scipy.optimize import fsolve
from dispRel_EBW import disp_rel as disp_rel_ebw
import derivForHamEq_EBW2 as ebw

Num= 1

def Z(mu):
   return (1 + 1j*mu)**-1

def get_ray(init_cond,ncycle):

    #init_cond=np.array([z0,x0,Nz0,Nx0])   
    def f(q,t):
        z,x,ReNz,ReNx,ImNz,ImNz = q.reshape(6,Num)
        # all the variables are in global CS
        r = np.array([z,x])
        N=np.array([ReNz,ReNx])
	
        X=p.X(r)*Z(1e-4)**2
        Y=mf.Yabs(r)*Z(1e-4)
        gamma = p.gamma(10)
        
        
	"""
        if localNx<1:
            vg = -np.real( dD_dN( X,Y,gamma,localNz,localNx,dNz_dN,dNx_dN  ) / dD_dw( X,Y,gamma,localNz,localNx  ) )

            a = np.real( dD_dr(X,Y,gamma,localNz,localNx,dNz_dr,dNx_dr,dX_dr,dY_dr,dgamma_dr) / dD_dw( X,Y,gamma,localNz,localNx  ) )
            A,B,C=get_ABC(X,Y,gamma,localNz,localNx)
            #print A*localNx**4
            arr = np.array([vg,a])
	"""

        dnz_dN=ri.dnz_dN(r,N)
        dnx_dN=ri.dnx_dN(r,N)
        dnz_dr=ri.dnz_dr(r,N)
        dnx_dr=ri.dnx_dr(r,N)
        dX_dr = p.dX_dr(r)
        dY_dr = mf.dY_dr(r)
        dgamma_dr = p.dgamma_dr(r)
        nz,nx=ri.nz(r,N),ri.nx(r,N)
      
        n = np.array([nz,nx])


        if nx<1:
            vg = -np.real( dD_dN( X,Y,gamma,nz,nx,dnz_dN,dnx_dN  ) / dD_dw( X,Y,gamma,nz,nx  ) )

            a = np.real( dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / dD_dw( X,Y,gamma,nz,nx  ) )
            Ima = np.imag( dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / dD_dw( X,Y,gamma,nz,nx  ) )
            print dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / dD_dw( X,Y,gamma,nz,nx  )  
            arr = np.array([vg,a,Ima])
            if conf.verbose ==True:
	    	print "dr",disp_rel(X,Y,gamma,nz,nx), x
        else:
            vg = -np.real( ebw.dD_dN( X,Y,gamma,nz,nx,dnz_dN,dnx_dN  ) / ebw.dD_dw( X,Y,gamma,nz,nx  ) )
       
            if vg[0]>1 or vg[1]>1:
                print "vg",vg[0],vg[1]

	    if conf.verbose==True:
            	print "EBW", x,X,disp_rel_ebw(X,Y,gamma,nz,nx), nx, nz
          
            a = np.real( ebw.dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / ebw.dD_dw( X,Y,gamma,nz,nx  ) )
            Ima= np.imag( ebw.dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / ebw.dD_dw( X,Y,gamma,nz,nx  ) )
	
	    arr = np.array([vg,a,Ima])
        return arr.reshape(6)


    def solv(init_cond,t):
        # initial cond
    
      
        init = init_cond
        print np.reshape(init,6*Num)
        soln = odeint(f,np.reshape(init,6*Num),t,full_output=0, printmessg=1) 


        return soln
    
    nstep= 20000 #number of steps is important. Too much can be bad apparently due too floating ariphmetics. 2e4 didnt work, 20 work well???

    

    t  = np.linspace(0, ncycle, nstep)
    soln=solv(init_cond,t)	
 
    return soln



if __name__ == '__main__':
	 
	z0 = np.array([0])
	x0 = np.array([1.01])#np.array([1.01])
	r = np.array([z0,x0])
	


	Nz0 = np.array([0.1])#np.array([ (1/np.sqrt(2) / (1 + 1/np.sqrt(2)))**0.5 - 0.1 ] )

	ncycle = 400
	params.Y=.8

	Y=mf.Yabs(r)
	X = p.X(r)
	Nx0 = coldNx_Xm(X,Y,Nz0)
	N = np.array([Nz0,Nx0])
	nz,nx=ri.nz(r,N),ri.nx(r,N)
	
	init_cond=np.array([z0,x0,Nz0,Nx0,np.array([0]),np.array([0])])  
	soln=get_ray(init_cond,ncycle)	
