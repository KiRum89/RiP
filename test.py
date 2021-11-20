import numpy as np
import matplotlib.pyplot as plt
from conf import p,mf
from multiprocessing import Pool
import sys
from dispRel import A,B,C, disp_rel
from coldDispRel import cold_disp_rel,coldNx_Xm,coldDispRel_pol,S,coldNx_Om,C
import refractive_index as ri
from scipy.optimize import brentq

import os
import hot as hot
import cold as cold
import params
from scipy.optimize import fsolve


Nz0 = np.array([0.1])#np.array([ (1/np.sqrt(2) / (1 + 1/np.sqrt(2)))**0.5 - 0.1 ] )
def hot_disp_rel_wrapp(nx,*data):
	X,Y,gamma,nz = data
	return np.real(disp_rel(X,Y,gamma,nz,nx))


ncycle =params.ncycle
#type of simulation


####
#f for pool

Te = params.Te #set the the Te

z0 = np.array([0])
x0=np.array([1.2])


r = np.array([z0,x0])
X = p.X(r)
Y = mf.Yabs(r)

data = X,Y,p.gamma(10),Nz0

Nx0 = coldNx_Xm(X,Y,Nz0)
Nx0=fsolve(hot_disp_rel_wrapp,Nx0,args=data)
print(Nx0)
N = np.array([Nz0,Nx0])
nz,nx=ri.nz(r,N),ri.nx(r,N)

init_cond=np.array([z0,x0,Nz0,Nx0])  
soln=hot.get_ray(init_cond,ncycle)	
#make gamma in the config!
print(init_cond)
soln = {"init_cond":init_cond,"soln":soln,"ncycle":ncycle,"gamma":p.gamma(10)}
                    
