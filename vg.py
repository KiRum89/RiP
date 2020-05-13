import refractive_index as ri
import conf as conf
from derivForHamEq import *
from hotDiTen import HotDiTen
from dispRel import A,B,C, disp_rel
from  scipy.integrate import odeint
from coldDispRel import cold_disp_rel,coldNx_Xm,coldDispRel_pol,S,coldNx_Om
from scipy.optimize import fsolve
from dispRel_EBW import disp_rel as disp_rel_ebw
import derivForHamEq_EBW2 as ebw
import sys


sys.path.insert(0,"/home/rumiantcev/RiP_general/plotting/")
from plot_Nx import start_Xm
## to check hot O mode in slab geometry


if conf.magnetic_field_model == "slab":
    print "magnetic slab"
    mf = __import__("mf_slab")

if conf.magnetic_field_model == "mirror":
    mf = __import__("mf_mirror")

if conf.magnetic_field_model == "simple":
    mf = __import__("mf_mirrorSimple")



if conf.plasma_model == "slab":
    print "plasma slab"
    p = __import__("plasma_slab")

if conf.plasma_model == "plasma":
    print "plasma"
    p = __import__("plasma")

doXm=True

Num= 1

gamma = p.gamma(10)    
def get_vg(r,N):
	X = p.X(r)
	Y = mf.Yabs(r)
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
            arr = np.array([vg,a])
        else:
            vg = -np.real( ebw.dD_dN( X,Y,gamma,nz,nx,dnz_dN,dnx_dN  ) / ebw.dD_dw( X,Y,gamma,nz,nx  ) )
	return vg

"""
soln = np.load("/home/rumiantcev/results/soln_2dXm.npy")
vg = np.zeros((soln[:,0].size,2))
for i in range(0,soln[:,0].size,10):
	j = i/10.	
	r =np.array([np.array([soln[j,0]]),np.array([soln[j,1]])])
	N = np.array([np.array([soln[j,2]]),np.array([soln[j,3]])]) 
	vg[j,:]=	get_vg(r,N).reshape(1,2)
"""
