import sys
import os
##set paths to packages
project = "RiP_general_new"
sys.path.insert(0,os.path.abspath(""))
import conf
from  scipy.integrate import odeint

import corrector



import refractive_index as ri
from conf import p,mf
import params

from derivForHamEq import *
from dispRel import A,B,C, disp_rel
from coldDispRel import cold_disp_rel,coldNx_Xm,coldDispRel_pol,S,coldNx_Om
from scipy.optimize import fsolve
from dispRel_EBW import disp_rel as disp_rel_ebw
import derivForHamEq_EBW2 as ebw

Num= 1

useCorrector = False; #some idea: use fsolve to find Nx from the dispersion relation if D drifts too far. Didnt work
Te=params.Te
def get_ray(init_cond,ncycle):

    #init_cond=np.array([z0,x0,Nz0,Nx0])   
    def f(q,t):
        z,x,Nz,Nx = q.reshape(4,Num)
        # all the variables are in global CS
        r = np.array([z,x])
        X=p.X(r)
        Y=mf.Yabs(r)
        gamma = p.gamma(Te)
	if useCorrector==True:

		data = X,Y,p.gamma(Te),Nz
		if np.abs(corrector.hot_disp_rel_wrapp(Nx,*data))>1e-4:
			sign = np.sign(np.abs(Nx))
			Nx=sign*fsolve(corrector.hot_disp_rel_wrapp,Nx,args=data)
	        
        N=np.array([Nz,Nx])
        

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
            vg = - dD_dN( X,Y,gamma,nz,nx,dnz_dN,dnx_dN  ) / dD_dw( X,Y,gamma,nz,nx  ) 

            a =  dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / dD_dw( X,Y,gamma,nz,nx  ) 
            if conf.verbose ==True:
	    	print('disp rel cold:{}, disp rel:{},x:{} '.format(cold_disp_rel(X,Y,nz,nx,"Xm"),disp_rel(X,Y,gamma,nz,nx), x))
            #if disp_rel(X,Y,gamma,nz,nx)>0:
	    #break
        else:
            vg = -np.real( ebw.dD_dN( X,Y,gamma,nz,nx,dnz_dN,dnx_dN  ) / ebw.dD_dw( X,Y,gamma,nz,nx  ) )
       
            if vg[0]>1 or vg[1]>1:
                print "vg",vg[0],vg[1]

	    if conf.verbose==True:
            	print( "EBW",'x:{},X:{},dispRel:{},nx:{},nz:{},'.format(x,X,disp_rel_ebw(X,Y,gamma,nz,nx), nx, nz))
	    #if  np.abs(disp_rel_ebw(X,Y,gamma,nz,nx))>1e-4:
	    #break         
            a = ebw.dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / ebw.dD_dw( X,Y,gamma,nz,nx  ) 

        arr = np.array([vg,a])
        arr = np.real(np.real(arr.astype(np.complex64)))    
        
        return arr.reshape(4*Num)


    def solv(init_cond,t):
        # initial cond
    
      
        init = np.real(init_cond)
        soln = odeint(f,np.reshape(init,4*Num),t,full_output=0, printmessg=1) 



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


	Y=mf.Yabs(r)
	X = p.X(r)
	Nx0 = coldNx_Xm(X,Y,Nz0)
	N = np.array([Nz0,Nx0])
	nz,nx=ri.nz(r,N),ri.nx(r,N)
	
	init_cond=np.array([z0,x0,Nz0,Nx0])  
	soln=get_ray(init_cond,ncycle)	
