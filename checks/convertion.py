import numpy as np
import matplotlib.pyplot as plt
import sys
import conf
"""
check what n_paral the SX-rays have at UHR after traversing the plasma in the wavechannel
"""
#soln1,soln1 set of solution with n_parallel in the oposite directions, func. "check" checks if the solution is physical ( the speed of light is not exceeded)

def nearUHR(soln,i):
	x_uhr = np.abs(soln[i,0,1]) #abs becasue it is in the opposite side
	#estimation of time when a ray reacehs the uhr

	print x_uhr
	return np.where(np.abs(soln[i,:,1]-x_uhr)<.01)[0][0]


