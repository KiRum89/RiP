import sys
import os
##set paths to packages
project = "RiP_general_new"
sys.path.insert(0,os.path.abspath(""))
import conf
from  scipy.integrate import odeint





import refractive_index as ri
from conf import p,mf


from derivForHamEq import *
from dispRel import A,B,C, disp_rel
from coldDispRel import cold_disp_rel,coldNx_Xm,coldDispRel_pol,S,coldNx_Om
from scipy.optimize import fsolve
from dispRel_EBW import disp_rel as disp_rel_ebw
import derivForHamEq_EBW2 as ebw




def f(r,N):
	
	X=p.X(r)
	Y=mf.Yabs(r)
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

	    a = np.real( dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / dD_dw( X,Y,gamma,nz,nx  ) )
	else:


	  
	    a = np.real( ebw.dD_dr(X,Y,gamma,nz,nx,dnz_dr,dnx_dr,dX_dr,dY_dr,dgamma_dr) / ebw.dD_dw( X,Y,gamma,nz,nx  ) )
	return a

