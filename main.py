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

z0 = np.array([0])
x0 = np.array([1.01])#np.array([1.01])
r = np.array([z0,x0])

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
doHarmonics=False
doSimpleMirror=True
doMirror=False

doOmodeBeam=False
doHarmonics_NzScan=False
doSXinBVP =False 

####
#f for pool

Te = params.Te #set the the Te
if doHarmonics==True:
	def f(Y):
		params.Y=Y
		X = p.X(r)
		Nz0 = np.array([np.sqrt(Y/(1+Y))])

		Nx0 = -coldNx_Om(X,Y,Nz0)#Om should be called "+" mode. we actually above the O-cutoff
		data = X,Y,p.gamma(Te),Nz0

		Nx0=fsolve(hot_disp_rel_wrapp,Nx0,args=data)

		N = np.array([Nz0,Nx0])
		nz,nx=ri.nz(r,N),ri.nx(r,N)

		init_cond=np.array([z0,x0,Nz0,Nx0])  
		
		print mf.Yabs(np.array([z0,x0])),init_cond
		soln=hot.get_ray(init_cond,ncycle)	

		soln = {"init_cond":init_cond,"soln":soln,"Y":params.Y,"ncycle":ncycle}
		#np.save(os.path.abspath("results/slab/Obeam/harm1/gamma10/soln_"+str(params.Y)),soln)
elif doSimpleMirror==True:
	def f(z0):
		z0 = np.array([z0])
		Nz0 = -np.array([0.03])
		#x0=-np.array([0.13])
		x0 = -np.array([0.20])

		
		r = np.array([z0,x0])
		X = p.X(r)
		Y = mf.Yabs(r)
		Nx0 = coldNx_Xm(X,Y,Nz0)
		data = X,Y,p.gamma(Te),Nz0
		
		Nx0=fsolve(hot_disp_rel_wrapp,Nx0,args=data)

		N = np.array([Nz0,Nx0])
		nz,nx=ri.nz(r,N),ri.nx(r,N)
		
		init_cond=np.array([z0,x0,Nz0,Nx0])  
		
		soln=hot.get_ray(init_cond,ncycle)	
		#make gamma in the config!
		soln = {"init_cond":init_cond,"soln":soln,"ncycle":ncycle,"gamma":p.gamma(Te)}
		#make file exist check
		np.save(os.path.abspath("results/mirror/simple/gamma"+str(Te) + "/soln_"+str(z0[0])),soln)
		return soln
elif doMirror==True:
	def f(z0):
		z0 = np.array([z0])
		Nz0 = np.array([-0.03])
		x0=-np.array([0.13])
		r = np.array([z0,x0])
		X = p.X(r)
		Y = mf.Yabs(r)
		Nx0 = coldNx_Xm(X,Y,Nz0)
		data = X,Y,p.gamma(10),Nz0
		
		Nx0=fsolve(hot_disp_rel_wrapp,Nx0,args=data)

		N = np.array([Nz0,Nx0])
		nz,nx=ri.nz(r,N),ri.nx(r,N)
		
		init_cond=np.array([z0,x0,Nz0,Nx0])  
		
		soln=hot.get_ray(init_cond,ncycle)	
		#make gamma in the config!
		soln = {"init_cond":init_cond,"soln":soln,"ncycle":ncycle,"gamma":p.gamma(10)}
		#make file exist check
		#np.save(os.path.abspath("results/mirror/simple/gamma100/soln_"+str(z0[0])),soln)
		return soln

elif doOmodeBeam==True:
	def f(z0):
		z0 = np.array([z0])
		x0=np.array([0.01])
		r = np.array([z0,x0])
		X = p.X(r)
		Y = mf.Yabs(r)
		Nz0 = np.sqrt(Y/(1+Y))-0.001
		print "Nz_opt",Nz0
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
		soln = {"init_cond":init_cond,"soln":soln,"ncycle":ncycle,"gamma":p.gamma(10)}
		#make file exist check

		#np.save(os.path.abspath("results/mirror/simple/Obeam/plasma_bvp/gamma10/soln_"+str(z0[0])),soln)
				
		return soln

elif doSXinBVP==True:
	import SXinFLiPS
	print "imported SXinFLiPS"
	
	r = np.array([0,-0.134]).reshape(2,1)
	def f(Nz0):
		return SXinFLiPS.f(r,Nz0)	
elif doHarmonics_NzScan==True:
	def f(Y):
		params.Y=Y
		Nz0_opt = np.array([np.sqrt(Y/(1+Y))])
		allSolns = []
		for Nz0 in np.linspace(Nz0_opt+0.001,Nz0_opt+0.001,1):
			Nz0 = np.array([Nz0])
			#TODO: use the cold dispersion relation to find which refractive index to use, becasue I want to do the angular scan. We need to choose the proper branch
			if (Nz0<Nz0_opt):
				#get the starting point of a ray				
				x0=np.array([brentq(cold_C_wrapp,1.001,40,args=(Y,Nz0))])+0.001
				print "x0",x0
			else: 
				x0=np.array([1.001])	#more or less arbitrary	
			z0 = np.array([0])		
			r = np.array([z0,x0])
			X = p.X(r)
			print r,X

			Nx0 = -coldNx_Om(X,Y,Nz0)
			print Nx0
			data = X,Y,p.gamma(10),Nz0
		
			Nx0=fsolve(hot_disp_rel_wrapp,Nx0,args=data)

			N = np.array([Nz0,Nx0])
			nz,nx=ri.nz(r,N),ri.nx(r,N)

			init_cond=np.array([z0,x0,Nz0,Nx0])  
			
			print mf.Yabs(np.array([z0,x0])),init_cond
			soln=hot.get_ray(init_cond,ncycle)	

			soln = {"init_cond":init_cond,"soln":soln,"Y":params.Y,"ncycle":ncycle}
			allSolns.append(soln)	
			#Idea: investigate vg_perp/v_parallel of as function of Nz in the slabb model. We willl find the that the dependy is weak! (hopefully) -> EBW is not affected much the injection angle		
		#np.save(os.path.abspath("results/slab/Obeam/harm1/gamma10/NzScan/soln_"+str(params.Y)),allSolns)

		return np.array(allSolns)
doPool =True 

if doPool:
        """	
	Y = mf.Yabs(r)

	Nzopt =np.sqrt(Y/(1+Y)) 
	"""
	p = Pool(40)
	##p.map(f,np.linspace(Nzopt+0.001,Nzopt+0.1,50))
	p.map(f,np.linspace(-0.1,0.1))

else:
	#Nz0=0.6687584399361054
	soln=f(0.8)#the input is not stanartiswd,can be Y or z0...
