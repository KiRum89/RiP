from __future__ import division
import numpy as np
from conf import mf,p
from matplotlib.pyplot import cm

import matplotlib.pyplot as plt

Num=500


Yz = np.zeros((500,500))
Yx = np.zeros((500,500))
X =  np.zeros((500,500))

i=0
z = np.linspace(-2,2,Num)
r = np.linspace(2,0,Num)
for x in r:
    x = x*np.ones(500)
    Yz[i,:] = mf.Yz((z,x))
    Yx[i,:] = mf.Yx((z,x))
    X[i,:] = p.X((z,x))
    i+=1


plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')

fig,(ax1)=plt.subplots(1,1)

fig.set_size_inches(5.4,5.5, forward=True)
#fig.set_figwidth(5.4,forward=True)
Y=np.sqrt(Yz**2+Yx**2)


im=ax1.imshow(X,cmap=plt.cm.Blues,extent=[-2,2,0,2] )


#cset=ax1.contour(z,r,Y,np.array([1/2,1]),linewidths=2, cmap = cm.Set2)
#ax1.clabel(cset,inline=True,fmt='Y=%1.f',fontsize=11)

cuhr=ax1.contour(z,r,np.sqrt(Y**2+X ),np.array([1]),linewidths=2, colors="red"  )
ax1.clabel(cuhr,inline=True,fmt="UHR",fontsize=12)

cO=ax1.contour(z,r,X ,np.array([1]),linewidths=2, colors="green"  )
ax1.clabel(cO,inline=True,fmt="O-cutoff",fontsize=12)


cR=ax1.contour(z,r,np.sqrt(0.25*Y**2+X )+0.5*Y,np.array([1]),linewidths=2, colors="black"  )
ax1.clabel(cR,inline=True,fmt='R',fontsize=12)

cL=ax1.contour(z,r,np.sqrt(0.25*Y**2+X )-0.5*Y,np.array([1]),linewidths=2, colors="black"  )
ax1.clabel(cL,inline=True,fmt='L-cutoff',fontsize=12)

ax1.streamplot(z,r,Yz,Yx,minlength=0.9, density=0.4,color="blue")    


ax1.set_xlabel(r"$z~\textrm{[m]} \parallel B_0 $",fontsize=12)
ax1.set_ylabel(r"$r~\textrm{[m]}$",fontsize=12)
ax1.set_aspect('equal')
#ax.set_xlim([-0.1,0.1])
#ax.set_ylim([-0.2,0.2])
#ax1.autoscale(False)

