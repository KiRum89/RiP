import numpy as np
import matplotlib.pyplot as plt
from local_N import Local_N
from plasma import Plasma


def checkNx(soln):

    p = Plasma(0.5,np.array([-1,1]))
    i = 0
    for x in soln[10345:10370,1]: 
        
        r = np.array([soln[i,0],soln[i,1]])
        
        
        locN=Local_N(p.Yvec(r),p.gradY(r), np.array([soln[i,2],soln[i,3]]) )
        print locN.Nx(), locN.Nz()
        plt.scatter(x,locN.Nx())
        plt.scatter(x,locN.Nz())
        
        i+=1
