from selector import H
from derivOfHotDiTen import DerivOfHotDiTen
from hotDiTen import HotDiTen
from dispRel import A,B,C



def dA_du(u,Nz,dielTen, derivOfDielTen):
    
    Kxx,Kyy,Kzz,Cxy,Cxz,Cyz = dielTen
    dKxx_du, dKyy_du, dKzz_du, dCxy_du, dCxz_du, dCyz_du = derivOfDielTen
    return 2*(Nz*dCxz_du+Cxz*H(u,"Nz"))+dKxx_du+2*Cxz*dCxz_du

def dB_du(u,Nz,dielTen, derivOfDielTen):
    Kxx,Kyy,Kzz,Cxy,Cxz,Cyz = dielTen
    dKxx_du, dKyy_du, dKzz_du, dCxy_du, dCxz_du, dCyz_du = derivOfDielTen

    
    T1 = Nz**3*dCxz_du+3*Nz**2*Cxz*H(u,"Nz")
    T2 = dKxx_du+dKzz_du+2*Cyz*dCyz_du+2*Cxz*dCxz_du #Nz**2
    T3 = Kxx+Kzz+Cyz**2+Cxz**2 # H
    T4 = Cxy*dCyz_du + Cyz*dCxy_du+Kyy*dCxz_du+Cxz*dKyy_du #Nz
    T5 = Cxy*Cyz + Kyy*Cxz # H

    
    T6 = -Kxx*(dKyy_du+dKzz_du) - (Kyy+Kzz)*dKxx_du + 2*Cxy*dCxy_du\
         -2*Kxx*Cyz*dCyz_du - Cyz**2*dKxx_du
    T7 = -2*(Cxz*Cyz*dCxy_du + Cxz*Cxy*dCyz_du + Cyz*Cxy*dCxz_du) - Cxz**2*dKyy_du - 2*Cxz*Kyy*dCxz_du
    #print "T",T1,T2,T3,T4,T5,T6,T7
    return 2*T1 + T2*Nz**2 + 2*Nz*T3*H(u,"Nz")  - 2*T4*Nz - 2*T5*H(u,"Nz")+T6 + T7


def dC_du(u,Nz,dielTen, derivOfDielTen):
    Kxx,Kyy,Kzz,Cxy,Cxz,Cyz = dielTen
    dKxx_du, dKyy_du, dKzz_du, dCxy_du, dCxz_du, dCyz_du = derivOfDielTen
    return C(Nz,dielTen)/Kzz*dKzz_du + Kzz*( 4*Nz**3*H(u,"Nz") - 2*Nz*(Kxx+Kyy)*H(u,"Nz") - Nz**2*(dKxx_du+dKyy_du)+Kxx*dKyy_du+Kyy*dKxx_du-2*Cxy*dCxy_du )



def dD_du(u,X,Y,gamma,Nz,Nx):



    num_term = 20 #in future will be global or define in a settings file
    
    hdt=HotDiTen(X,Y,gamma,Nz,Nx,num_term)
    dielTen = [hdt.Kxx(),hdt.Kyy(),hdt.Kzz(),hdt.Cxy(),hdt.Cxz(),hdt.Cyz()]
    ddt = DerivOfHotDiTen(X,Y,gamma,Nz,Nx,num_term)
    derivOfDielTen = [ddt.dKxx_du(u), ddt.dKyy_du(u), ddt.dKzz_du(u), ddt.dCxy_du(u), ddt.dCxz_du(u), ddt.dCyz_du(u)]
   

    #print dielTen
    #print dA_du(u,Nz,dielTen,derivOfDielTen),dB_du(u,Nz,dielTen,derivOfDielTen),dC_du(u,Nz,dielTen,derivOfDielTen)
    return Nx**4*dA_du(u,Nz,dielTen,derivOfDielTen) + 4*Nx**3*A(Nz,dielTen)*H(u,"Nx") + Nx**2*dB_du(u,Nz,dielTen,derivOfDielTen) + 2*Nx*B(Nz,dielTen)*H(u,"Nx") + dC_du(u,Nz,dielTen,derivOfDielTen)

 
