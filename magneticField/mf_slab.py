import numpy as np
import conf
import params


def q(r):
    return 1/Yabs(r)


def q_vec(r):
    return 1/np.array([Yz(r),Yx(r)])


def Yvec(r):
    return np.array([Yz(r),Yx(r)])

def Yabs(r):
    return np.sqrt(np.sum(Yvec(r)**2,axis=0))


def Yz(r):
    z,x = r
    return np.array([params.Y])#1./np.sqrt(2)*np.ones(z.size)

def Yx( r):
    z,x=r
    return np.zeros(z.size)


def gradY(r):
    z,x = r
    dr = 1e-6
    r_dx = (z,x+dr)
    r_dz = (z+dr,x)

    dY_dz=Yvec(r_dz) - Yvec(r)
    dY_dx=Yvec(r_dx) - Yvec(r)
    return np.array([dY_dz,dY_dx])/dr


def dY_dr(r):
    return np.sum(gradY(r)*Yvec(r),axis=1)/Yabs(r)

def grad_q(r):
    z,x = r
    dr = 1e-6
    r_dx = np.array([z,x+dr])
    r_dz = np.array([z+dr,x])
    dq_dz=q_vec(r_dz) - q_vec(r)
    dq_dx=q_vec(r_dx) - q_vec(r)
    return np.array([dq_dz,dq_dx])/dr

def q_abs(r):
    return np.sqrt( np.sum(q_vec(r)**2) )
