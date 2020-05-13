from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def R(X,Y):
    return 1 - X/(1 - Y)

def L(X,Y):
    return 1 - X/(1 + Y)

def S(X,Y):
    return 0.5*( R(X,Y)+L(X,Y) )

def B(X,Y,gamma):
    
    return X/Y**2*gamma**2 * ( 1 /(Y**2 - 1 ) - 1/(4*Y**2 - 1) )

def n1_perp2(X,Y,gamma):

    return ( S(X,Y) + np.sqrt( S(X,Y)**2 - 4*R(X,Y)*L(X,Y)*B(X,Y,gamma) ) )/2/B(X,Y,gamma)

def n2_perp2(X,Y,gamma):

    return ( S(X,Y) - np.sqrt( S(X,Y)**2 - 4*R(X,Y)*L(X,Y)*B(X,Y,gamma) ) )/2/B(X,Y,gamma)


def X(x):
    p1=2
    p2=3
    X1=0.001
    X0=2
    result = X0 + (X1-X0)*( 1-( np.abs(x)/xb )**p1 )**p2
    

    return result

xb = 4*np.pi


"""
ax = plt.subplot(111)
gamma = 0.006


X1 = np.linspace(0.01,4,100000)
X2 = np.linspace(0.01,4,10000)
Y=0.9

ax.plot(X1, n1_perp2(X1,Y,gamma), label="+")
ax.plot(X2, n2_perp2(X2,Y,gamma),label='-')

#ax.plot(X1, np.abs(np.imag(n1_perp2(X1,Y,gamma))), label="+")
ax.plot(X2, np.abs(np.imag(n2_perp2(X2,Y,gamma))),label='-')
#ax.set_yscale('log')

ax.legend()


X3 = np.linspace(0.1,0.8,100000)
ax.plot(X3, n2_perp2(X3,Y,gamma))

gamma = 0.062561
ax.plot(X1, n1_perp2(X1,Y,gamma))
ax.plot(X2, n2_perp2(X2,Y,gamma))
ax.plot(X3, n2_perp2(X3,Y,gamma))


X3 = np.linspace(0.01,2,100000)
gamma = 0.62561
ax.plot(X1, n1_perp2(X1,Y,gamma))
ax.plot(X2, n2_perp2(X2,Y,gamma))
ax.plot(X3, n2_perp2(X3,Y,gamma))


Xuhr=0.631
print "fdfsdf" ,S(Xuhr,Y)
gamma = 0.0062561
ax.scatter( Xuhr,np.sqrt( R(Xuhr,Y)*L(Xuhr,Y) / B(Xuhr,Y,gamma) ) )
plt.figure()

XX = np.linspace(0,4,10000)
plt.plot(XX,S(XX,Y)**2, label="S")
plt.plot(XX,4*R(XX,Y)*L(XX,Y)*B(XX,Y,gamma), label="RLB")
plt.legend()

#plt.yscale('log')

"""

def D(X,Y,gamma,n_perp):
    return B(X,Y,gamma)*n_perp**4 - S(X,Y)*n_perp**2 + R(X,Y)*L(X,Y)

def dD_dx(X,Y,gamma,n_perp):
    dx = 1e-6
    dn = 1e-6
    return (D(X+dx,Y,gamma,n_perp) - D(X,Y,gamma,n_perp))/dx * 1


def dD_dn(X,Y,gamma,n_perp):

    dn = 1e-6
    return (D(X,Y,gamma,n_perp+dn) - D(X,Y,gamma,n_perp))/dn 


def dD_dw(X,Y,gamma,n_perp):
    dx = 1e-6
    dX_dw = - 2*X
    dY_dw =  -Y

    return (D(X+dx,Y,gamma,n_perp) - D(X,Y,gamma,n_perp))/dx * dX_dw\
    +(D(X,Y+dx,gamma,n_perp) - D(X,Y,gamma,n_perp))/dx * dY_dw\
    -(D(X,Y,gamma,n_perp+dx) - D(X,Y,gamma,n_perp))/dx * n_perp


"""
fig = plt.figure()
ax1 = plt.subplot(111)
gamma = 0.0062561


X1 = np.linspace(0.64,4,1000000)
n_perp1=np.sqrt(n2_perp2(X1,Y,gamma))
vg1 = -dD_dn(X1,Y,gamma,n_perp1.reshape(1,X1.size))/dD_dw(X1,Y,gamma,n_perp1.reshape(1,n_perp1.size))*3e8

X2 = np.linspace(0.1,0.41,1000)
n_perp2=np.sqrt(n2_perp2(X2,Y,gamma)).reshape(1,X2.size)
vg2 = -dD_dn(X2,Y,gamma,n_perp2)/dD_dw(X2,Y,gamma,n_perp2)*3e8

X3 = np.linspace(0.64,0.8,1000)
n_perp3=-np.sqrt(n2_perp2(X3,Y,gamma)).reshape(1,X3.size)
vg3 = -dD_dn(X3,Y,gamma,n_perp3)/dD_dw(X3,Y,gamma,n_perp3)*3e8

ax1.plot(X1,vg1[0])
ax1.plot(X2,vg2[0])
ax1.plot(X3,vg3[0])
ax1.set_yscale('log')
"""

gamma = 0.0062561
Z = 1e-3
Y = 0.7


"""
fig = plt.figure()
ax1 = plt.subplot(111)
gamma = 0.0062561
Z = 1e-3
Y = 0.7

X = np.linspace(0.1,1.7,1000)
n1=n1_perp2(X,Y,gamma)



n2=n2_perp2(X,Y,gamma)

plt.plot(X,n1, label='n1')
plt.plot(X,n2,label='n2')
plt.yscale('log')
plt.legend()

plt.figure()
"""

xaxb= np.linspace(4,4*np.pi,5000)
X1 = X(xaxb)
n_perp1=-np.sqrt(n2_perp2(X1,Y,gamma))
vg1 = -dD_dn(X1,Y,gamma,n_perp1.reshape(1,X1.size))/dD_dw(X1,Y,gamma,n_perp1.reshape(1,n_perp1.size))
#plt.plot(X1,vg1.reshape(X1.size,1))


def phase_real(xaxb):
    arr = []
    for i in range(0,xaxb.size):
        
        arr.append(np.sum(n_perp1[0:i]*np.diff(xaxb)[0]))
    
    return np.array(arr)

def phase_imag(xaxb):
    arr = []
    for i in range(0,xaxb.size):
        
        arr.append(np.sum( -Z*( 1/vg1.reshape(1,X1.size)[0][0:i] - n_perp1[i])*np.diff(xaxb)[0]) )
    return np.array(arr)


plt.plot(xaxb,np.cos(phase_real(xaxb))*np.exp(phase_imag(xaxb)))

