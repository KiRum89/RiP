import numpy as np
from  hotDiTen import HotDiTen


def A(Nz,dielTen):
    Kxx,Kyy,Kzz,Cxy,Cxz,Cyz = dielTen
    return 2*Cxz*Nz+Kxx+Cxz**2

def B(Nz,dielTen):
    Kxx,Kyy,Kzz,Cxy,Cxz,Cyz = dielTen
    #print "B terms wich should be 0",2*Cxz, Cyz**2 + Cxz**2, Cxy*Cyz + Kyy*Cxz,Kxx*Cyz**2,2*Cxz*Cyz*Cxy,  Cxz**2*Kyy
    
    return 2*Cxz*Nz**3 + ( Kxx + Kzz + Cyz**2 + Cxz**2 )*Nz**2 \
        -2*( Cxy*Cyz + Kyy*Cxz )*Nz - Kxx*(Kyy+Kzz)\
        +Cxy**2 - Kxx*Cyz**2 - 2*Cxz*Cyz*Cxy - Cxz**2*Kyy

def C(Nz,dielTen):
    Kxx,Kyy,Kzz,Cxy,Cxz,Cyz = dielTen
    return Kzz*( Nz**4 - ( Kxx + Kyy ) * Nz**2 + Kxx*Kyy - Cxy**2 )

def disp_rel(X,Y,gamma,Nz,Nx): 


    num_term = 20 #in future will be global or define in a settings file
    
    hdt=HotDiTen(X,Y,gamma,Nz,Nx,num_term)
    dielTen = [hdt.Kxx(),hdt.Kyy(),hdt.Kzz(),hdt.Cxy(),hdt.Cxz(),hdt.Cyz()]

    return A(Nz,dielTen)*Nx**4 + B(Nz,dielTen)*Nx**2 + C(Nz,dielTen)

