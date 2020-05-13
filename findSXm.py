import numpy as np
from coldDispRel import coldDispRel_pol
from plasma import SlabPlasma
from coldDispRel import cold_disp_rel

    
def findSX(soln):

    xcutInd=np.where(np.abs(soln[:,1]>0.9))[0][-1]

    xcut = soln[xcutInd,1]
    zcut=soln[xcutInd,0]
    rcut = np.array([ zcut,1 ])

    Nz = soln[xcutInd,2] 
    Nx = -soln[xcutInd,3] 
    N = np.array([Nz,Nx])
    p = SlabPlasma(1)

    while (np.abs(cold_disp_rel(p.X([zcut,xcut]),p.Yabs([zcut,xcut]),Nz,Nx,"Xm"))>1e-5):
        print cold_disp_rel(p.X([zcut,xcut]),p.Yabs([zcut,xcut]),Nz,Nx,"Xm"), rcut[1]
        xcut +=1e-6

    return [xcut,N]
#[1.0258487913687231, array([ 0.6       ,  0.08655108])]
