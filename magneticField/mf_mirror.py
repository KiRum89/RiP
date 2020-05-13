import numpy as np
from scipy.special import ellipe
from scipy.special import ellipk
from scipy.constants import mu_0,e,m_e



a = 0.38
coils_out = np.array([-0.63,0.63])



def q(r):
    return 1/Yabs(r)


def q_vec(r):
    return 1/np.array([Yz(r),Yx(r)])


def Yvec(r):
    return np.array([Yz(r),Yx(r)])

def Yabs(r):
    return np.sqrt(np.sum(Yvec(r)**2,axis=0))

freq = 2.45e9

def Yz_out( r):
    z,x = r
    R = np.abs(x/ a)
    Z = (1/ a*( coils_out[:,np.newaxis]-z))
    k = ( 4*R/((1+R)**2 + Z**2) )**0.5
    T1 = 1/np.sqrt((1+R)**2 + Z**2)
    T2 = (1-R**2-Z**2)/( (1-R)**2+Z**2 )*ellipe(k**2)+ellipk(k**2)
    #return e/m_e*(mu_0*80*I/2/np.pi)/2/np.pi/freq/a*np.sum(T1*T2,axis=0)
    I1 = 500
    A1 = e/m_e*(mu_0*80*I1/2/np.pi)/2/np.pi/freq/a
    left_coil=A1*T1[0,:]*T2[0,:]

    I2 = 500
    A2 = e/m_e*(mu_0*80*I2/2/np.pi)/2/np.pi/freq/a
    right_coil=A2*T1[1,:]*T2[1,:]
    return left_coil + right_coil

def Yx_out( r):
    z,x = r

    R = np.abs(x/ a)
    Z = 1/ a*( coils_out[:,np.newaxis]-z)
    k = ( 4*R/((1+R)**2 + Z**2) )**0.5
    T1 = Z/R/np.sqrt((1+R)**2 + Z**2)
    T2 = (1+R**2+Z**2)/( (1-R)**2+Z**2 )*ellipe(k**2)-ellipk(k**2)

    I1 = 500
    A1 = -e/m_e*(mu_0*80*I1/2/np.pi)/2/np.pi/freq/a*np.sign(x)
    left_coil=A1*T1[0,:]*T2[0,:]

    I2 = 500
    A2 = -e/m_e*(mu_0*80*I2/2/np.pi)/2/np.pi/freq/a*np.sign(x)
    right_coil=A2*T1[1,:]*T2[1,:]
    return left_coil + right_coil


coils_in = np.array([-0.34,-0.12,0.12,0.34])

def Yz_in(r):
    z,x = r
    I = 300
    R = np.abs(x/ a)
    Z = (1/ a*( coils_in[:,np.newaxis]-z))
    k = ( 4*R/((1+R)**2 + Z**2) )**0.5
    T1 = 1/np.sqrt((1+R)**2 + Z**2)
    T2 = (1-R**2-Z**2)/( (1-R)**2+Z**2 )*ellipe(k**2)+ellipk(k**2)
    return e/m_e*(mu_0*40*I/2/np.pi)/2/np.pi/freq/a*np.sum(T1*T2,axis=0)

def Yx_in(r):
    z,x = r
    I = 300
    R = np.abs(x/ a)
    Z = 1/ a*( coils_in[:,np.newaxis]-z)
    k = ( 4*R/((1+R)**2 + Z**2) )**0.5
    T1 = Z/R/np.sqrt((1+R)**2 + Z**2)
    T2 = (1+R**2+Z**2)/( (1-R)**2+Z**2 )*ellipe(k**2)-ellipk(k**2)
    return -e/m_e*(mu_0*40*I/2/np.pi)/2/np.pi/freq/a*np.sum(T1*T2,axis=0)*np.sign(x)



def Yz(r):
    return ( Yz_out(r)+Yz_in(r))

def Yx(r):
    return ( Yx_out(r)+Yx_in(r))




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

def flux(r):
    z,x = r
    R = np.abs(x/a)
    Z = 1/a*(coils[:,np.newaxis]-z)
    k = ( 4*R/((1+R)**2 + Z**2) )**0.5
    T1 = R/np.sqrt((1+R)**2 + Z**2)
    T2 = ((2-k**2)*ellipk(k**2)-2*ellipe(k**2))/k**2
    return np.sum(T1*T2,axis=0)/1.5
