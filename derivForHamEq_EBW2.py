from selector import H
from derivOfHotDiTen import DerivOfHotDiTen
from hotDiTen import HotDiTen
from dispRel import A,B,C
from derivOfDispRel import dA_du, dB_du, dC_du
import numpy as np

def dD_du(u,X,Y,gamma,Nz,Nx):



    num_term = 20 #in future will be global or define in a settings file
    
    hdt=HotDiTen(X,Y,gamma,Nz,Nx,num_term)
    dielTen = [hdt.Kxx(),hdt.Kyy(),hdt.Kzz(),hdt.Cxy(),hdt.Cxz(),hdt.Cyz()]
    ddt = DerivOfHotDiTen(X,Y,gamma,Nz,Nx,num_term)
    derivOfDielTen = [ddt.dKxx_du(u), ddt.dKyy_du(u), ddt.dKzz_du(u), ddt.dCxy_du(u), ddt.dCxz_du(u), ddt.dCyz_du(u)] 


    val=dA_du(u,Nz,dielTen,derivOfDielTen) + dB_du(u,Nz,dielTen,derivOfDielTen)/Nx**2 - 2/Nx**3*B(Nz,dielTen)*H(u,"Nx") + dC_du(u,Nz,dielTen,derivOfDielTen)/Nx**4 - 4*C(Nz,dielTen)/Nx**5*H(u,"Nx")


    #print val,derivOfDielTen

    return val



def dD_dN(X,Y,gamma,Nz,Nx,dNz_dN,dNx_dN):
    val=dD_du("Nx",X,Y,gamma,Nz,Nx)*dNx_dN+dD_du("Nz",X,Y,gamma,Nz,Nx)*dNz_dN
    if np.abs(val)[0] < 1e-12 or np.abs(val)[1] < 1e-12:
        print( "dD_dN is small")
    return val



def dD_dr(X,Y,gamma,Nz,Nx,dNz_dr,dNx_dr,dX_dr,dY_dr,dgamma_dr):

    val=dD_du("X",X,Y,gamma,Nz,Nx)*dX_dr + dD_du("Y",X,Y,gamma,Nz,Nx)*dY_dr + dD_du("Nx",X,Y,gamma,Nz,Nx)*dNx_dr + dD_du("Nz",X,Y,gamma,Nz,Nx)*dNz_dr + dD_du("gamma",X,Y,gamma,Nz,Nx) * dgamma_dr

    #if   np.abs(val)[0] < 1e-12 or np.abs(val)[1] < 1e-12:
    #    print "dD_dr is small"
    
    return val

#removed 1/omega from the equation because it will cancel in the ray tracing equation
def dD_dw(X,Y,gamma,Nz,Nx):

    val=-( 2*X*dD_du("X",X,Y,gamma,Nz,Nx)+Y*dD_du("Y",X,Y,gamma,Nz,Nx)+Nx*dD_du("Nx",X,Y,gamma,Nz,Nx)+Nz*dD_du("Nz",X,Y,gamma,Nz,Nx) )

    #if np.abs(val)<1e-12:
    #    print "dD_dw is small"

    return val

