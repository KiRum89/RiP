import warnings
import numpy as np
import matplotlib.pyplot as plt
from derivForHamEq import *
from plasma import SlabPlasma, Plasma
from local_N import Local_N
vgx = []
vgz = []
lNz = []




warnings.filterwarnings('always')
warnings.simplefilter("always")

soln = np.load("/home/kirill/wendsday_seminar/trajectory_demo_trapping.npy")


def group_velocity(soln):

    
    coils = np.array([-0.63,0.63])
    p = Plasma(0.3,coils)
        


    
                 
    for i in range(0,582):
        #soln[:,1].size,10):


        print i
        z = soln[i,0]
        x = soln[i,1]
        #p = Plasma(1)
        X = p.X([z,x])
    
        Y = p.Yabs([z,x])
        gamma = p.gamma(1)
        Nz = soln[i,2]
        Nx = soln[i,3]



        localN = Local_N(p.Yvec([z,x]),p.gradY([z,x]),np.array([Nz,Nx]))

        dNz_dN = localN.dNz_dN()
        dNx_dN=localN.dNx_dN()     
        localNz=localN.Nz()
        localNx = localN.Nx()
        vgz.append((dD_dN( X,Y,gamma,localNz,localNx,dNz_dN,dNx_dN   )[0] / dD_dw(  X,Y,gamma,localNz,localNx   ) ) )
        vgx.append((dD_dN( X,Y,gamma,localNz,localNx,dNz_dN,dNx_dN   )[1] / dD_dw( X,Y,gamma,localNz,localNx  ) ) )
        lNz.append(localN.Nz())
        #print z,x,X,Y,Nz,Nx#dD_dN( X,Y,gamma,localNz,localNx,dNz_dN,dNx_dN   )[0] , dD_dw(  X,Y,gamma,localNz,localNx   ) 
    return np.array(vgz),np.array(vgx)





vg = group_velocity(soln)
