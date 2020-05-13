import numpy as np
import sys




import conf as conf
import matplotlib.pyplot as plt
import derivForHamEq_c as der
from coldDispRel import cold_disp_rel, cold_disp_rel_Nx2, coldNx_Xm, coldNx_Om
import refractive_index as ri
import sys
sys.path.insert(0,"/home/rumiantcev/RiP_general/plotting")
import utils as u


if conf.magnetic_field_model == "slab":
    mf = __import__("mf_slab")

if conf.magnetic_field_model == "mirror":
    mf = __import__("mf_mirror")

if conf.magnetic_field_model == "simple":
    mf = __import__("mf_mirrorSimple")


if conf.plasma_model == "slab":
    p = __import__("plasma_slab")

if conf.plasma_model == "plasma":
    print "plasma"
    p = __import__("plasma")


#Y,Nz should be arrays 1x1
def start_Xm(Y,Nz):
	
	
	X = np.linspace(1.0001,3,1000)

	 
	Nx0_Xm=coldNx_Xm(X,Y,Nz)
	Nx0_Om=coldNx_Om(X,Y,Nz)

	arrNan=np.where(np.isnan(Nx0_Xm)==True)[0]
	Nx_turnP = coldNx_Om(X[arrNan[0]-1],Y,Nz)		 
	ind1= np.where(Nx0_Om<Nx_turnP)[0]
	
	plus_mode_ind=ind1	
	if plus_mode_ind.size!=0:
		return X[plus_mode_ind[0]]

	else:
		return X[arrNan[0]] 
			
		
plot=True
if plot==True:

	z0 = np.zeros(1000)
	x0 = np.linspace(0,2,1000)



	r = np.array([z0,x0])


	Y=mf.Yabs(r)
	X = p.X(r)



	import params

	Y = params.Y*np.ones(1000)
	Nz0 = u.Nz_opt(Y[0])+0.1

	Nx0_Xm=coldNx_Xm(X,Y,Nz0)
	Nx0_Om=coldNx_Om(X,Y,Nz0)


	plt.plot(X,Nx0_Om)
	plt.plot(X,Nx0_Xm)
	plt.scatter(start_Xm(Y[0],Nz0),0)


