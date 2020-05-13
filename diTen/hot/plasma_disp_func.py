
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
from scipy.special import wofz


def Z(x):
    return wofz(x)*1j*np.pi**0.5

## first derivative of wofz (analytically)
def Zp(x):
    return -2 - 2*x*Z(x)

def dawsn_expansion(x):
    # Accurate to order x^-9, or, relative to the first term x^-8
    # So when x > 100, this will be as accurate as you can get with
    # double floating point precision.
    y = 0.5 * x**-2
    return 1/(2*x) * (1 + y * (1 + 3*y * (1 + 5*y * (1 + 7*y))))

def dawsn_expansion_drop_first(x):
    y = 0.5 * x**-2
    return 1/(2*x) * (0 + y * (1 + 3*y * (1 + 5*y * (1 + 7*y))))

def dawsn_expansion_drop_first_two(x):
    y = 0.5 * x**-2
    return 1/(2*x) * (0 + y * (0 + 3*y * (1 + 5*y * (1 + 7*y))))

def blend(x, a, b):
    # Smoothly blend x from 0 at a to 1 at b
    y = (x - a) / (b - a)
    y *= (y > 0)
    y = y * (y <= 1) + 1 * (y > 1)
    #print y*y * (3 - 2*y)
    return y*y * (3 - 2*y) # why we cant simply write y?

def g(x):
    """Calculate `x + (1-2x^2) D(x)`, where D(x) is the dawson function"""
    # For x < 50, use dawsn from scipy
    # For x > 100, use dawsn expansion
    b = blend(x, 50, 100)
    y1 = x + (1 - 2*x**2) * special.dawsn(x)
    y2 = dawsn_expansion_drop_first(x) - dawsn_expansion_drop_first_two(x) * 2*x**2
    return b*y2 + (1-b)*y1

def Zpp(x):
    # only return the imaginary component
    return -2j*np.pi**0.5 * g(x)


"""
x = np.logspace(0, 5, 2000)
dx = 1e-3
plt.plot(x, (Zp(x+dx) - Zp(x-dx)).imag/(2*dx))
plt.plot(x, Zpp(x).imag)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
"""
