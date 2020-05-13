import numpy as np
import conf as conf

from conf import mf,p



def nz(r,N):
    
    return np.sum( N*mf.Yvec(r)/mf.Yabs(r),axis=0 )

def nx(r,N):
    return np.sqrt(np.sum(N**2,axis=0)-nz(r,N)**2)


def dnz_dN(r,N):
    return mf.Yvec(r)/mf.Yabs(r)

def dnx_dN(r,N):
    return (N-nz(r,N)*mf.Yvec(r)/mf.Yabs(r))/nx(r,N)

    
def dnz_dr(r,N):
    return np.sum(mf.gradY(r)*nx(r,N)/mf.Yabs(r)*dnx_dN(r,N),axis=1 ) #np.array([0,0])

def dnx_dr(r,N):
        return -nz(r,N)/nx(r,N)*dnz_dr(r,N)
