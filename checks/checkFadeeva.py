import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz,iv, ivp


def Z(x):
    return 1j*np.pi**0.5*wofz(x)

def Zp(x):
    return -2 - 2*x*Z(x)

def Zpp(x):
    return (Z(x)+x*Zp(x))*x
x = np.float64(np.linspace(1e4,14e4,1000))
#plt.plot(x, Zpp(x),"o")
plt.plot(x[:-1],np.diff(Zp(x))/np.diff(x))
#plt.plot(x,x*Zp(x))

#plt.plot(x,Z(x))
