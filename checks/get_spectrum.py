import numpy as np
from isPhysical import check
"""
idea: what if we get k-spectrum at O-cutoff of the secondary SX-wave. Or at UHR. How diferent will it be from the original k-spectrum?
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.special import zeta
import numpy as np
import os
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import sys
import os
sys.path.insert(0,os.path.abspath(""))
sys.path.insert(0,os.path.abspath("plotting"))

from conf import mf,p
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import style
pathGamma = "gamma10"
path = os.path.abspath("results/mirror/simple/Obeam/gamma10")
#
from flips import get_part
def isSXwave(nx):
	if np.abs( nx)<5:
		return True
	else:
		return False
soln=np.load(path+"/soln2.npy")

soln = np.array(soln.item()["soln"])
X,Y,z,r = get_part((-0.6,0.6),(-0.25,0.25))
ind_uhr=np.where(np.abs(X[:,250]+Y[:,250]**2-1)<0.02)[0][0]
d=[]
r_uhr = np.abs(r[ind_uhr])-0.01
def spectrum(soln):
	ind2=[]
	for i in range (0,50):
		if check(soln,i):
			ind=np.where(np.abs(soln[i,:-5000,3])<1)[0]
			d.append( np.min(np.abs(soln[i,100:-5000,1]-r_uhr)))
			
			#ind2.append(soln[i,,3])
			print soln[i,np.argmin(np.abs(soln[i,100:,1]+100-r_uhr)),2]
			ind2.append(soln[i,np.argmin(np.abs(soln[i,ind,1]-r_uhr)),1])

	return np.array(ind2)			
solnz=[]
solnx=[]
t=10000
def spectrum2(soln):
	soln = soln[:,0:t,:]
	ind2=[]
	for i in range (0,50):
		if check(soln,i):
			try:
				print i
				soln2 = soln[i,soln[i,:,1]>0,:]
				print soln2
				solnz.append(soln2[np.argmin(np.abs(soln2[:,3])),2])
				solnx.append(soln2[np.argmin(np.abs(soln2[:,3])),3])


			except IndexError:
				print 12

	return solnz,solnx 
if __name__ == '__main__':
	soln=soln[np.argsort(soln[:,0,0]),:,:] 
	spect = spectrum2(soln)
