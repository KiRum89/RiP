import numpy as np
from dispRel import A,B,C, disp_rel

def hot_disp_rel_wrapp(nx,*data):
	X,Y,gamma,nz = data
	return np.real(disp_rel(X,Y,gamma,nz,nx))

