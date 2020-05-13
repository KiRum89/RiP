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
#when C=0, is the genereliside reflection condtion. We will use to find the start of SX-wave to investigate O-X-B convretion
def cold_C_wrapp(X,*data):
	Y,nz = data
	return C(X,Y,nz)


ncycle =params.ncycle
#type of simulation
doSimpleMirror=True


####
#f for pool

Te = params.Te #set the the Te
doOmodeBeam=True

if doOmodeBeam==True:
	def f(delta_Nz0):

                z0 = np.array([0])
		x0=np.array([1.05])

		r = np.array([z0,x0])
		X = p.X(r)
		Y = mf.Yabs(r)
		Nz0 = np.sqrt(Y/(1+Y))+delta_Nz0
		print "Nz_opt",X,Y
		Nx0 = coldNx_Om(X,Y,Nz0)
                
		data = X,Y,p.gamma(10),Nz0
		Nx0=fsolve(hot_disp_rel_wrapp,Nx0,args=data)
		print Nx0
		N = np.array([Nz0,Nx0])
		nz,nx=ri.nz(r,N),ri.nx(r,N)
		
		init_cond=np.array([z0,x0,Nz0,Nx0])  
		print init_cond	
		soln=hot.get_ray(init_cond,ncycle)	
		#make gamma in the config!
                print(init_cond)
		soln = {"init_cond":init_cond,"soln":soln,"ncycle":ncycle,"gamma":p.gamma(10)}
		#make file exist check

				
		return soln


doPool = False 


if doPool:
	p = Pool(40)
	p.map(f,np.linspace(-0.01,0.01))
else:
	soln=f(0.01)#the input is not stanartiswd,can be Y or z0...
