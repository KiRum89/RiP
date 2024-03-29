import numpy as np
from scipy.special import ellipe
from scipy.special import ellipk






def q(r):
    return 1/Yabs(r)


def q_vec(r):
    return 1/np.array([Yz(r),Yx(r)])


def Yvec(r):
    return np.array([Yz(r),Yx(r)])

def Yabs(r):
    return np.sqrt(np.sum(Yvec(r)**2,axis=0))


def Yz( r):
    z,x = r

    #return np.array([(z/0.3612)**2+0.85]) #parabolic fir for mf_ebw_only
    return (z/0.3)**2+0.7*np.ones(z.size)#(z/4)**2+0.7 #6coils simple

def Yx( r):
    z,x=r
    return np.zeros(z.size)# check first without Yx-np.sum(T1*T2,axis=0)*np.sign(x)

def gradY(r):
    z,x = r
    dr = 1e-6
    r_dx = (z,x+dr)
    r_dz = (z+dr,x)

    dY_dz=Yvec(r_dz) - Yvec(r)
    dY_dx=Yvec(r_dx) - Yvec(r)
    #return np.array([dY_dz,dY_dx])/dr
    return np.array([ np.array([2/0.3*z,np.zeros(z.size)]),np.array([np.zeros(z.size) ,np.zeros(z.size) ])])




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



