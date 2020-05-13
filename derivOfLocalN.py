

def dNz_dN(Yvec,Yabs):
    return Yvec/Yabs

def dNx_dN(Yvec,Yabs,locNz,locNx,globN):
    return ( globN - Yvec/Yabs*Nz )/Nx

def dNz_dr(q):
    return np.dot( gradY(q),Nx(q)/Yabs(q)*dNx_dN(q) ) #np.array([0,0])

def dNx_dr(q):
    #print -Nz(q)/Nx(q)*dNz_dr(q)
    return -Nz(q)/Nx(q)*dNz_dr(q)#

