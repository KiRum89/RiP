import numpy as np
def check(soln,i):



	vg=np.sqrt(np.diff(soln[i,:,0])**2 +np.diff(soln[i,:,1])**2 )
	dt = 600./2e4
	arr=np.where(np.abs(vg)/dt>1)
	if arr[0].size>0:
		return False	
		
	return True 

