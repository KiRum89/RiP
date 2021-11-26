from __future__ import division
import numpy as np
from conf import mf,p
from matplotlib.pyplot import cm

import matplotlib.pyplot as plt

def plot(soln):
    Num=500


    Yz = np.zeros((500,500))
    Yx = np.zeros((500,500))
    X =  np.zeros((500,500))

    i=0
    z = np.linspace(-2,2,Num)
    r = np.linspace(-1,1,Num)
    for x in r:
        x = x*np.ones(500)
        Yz[i,:] = mf.Yz((z,x))
        Yx[i,:] = mf.Yx((z,x))
        X[i,:] = p.X((z,x))
        i+=1



    fig,(ax1)=plt.subplots(1,1)

    fig.set_size_inches(5.4,5.5, forward=True)
    #fig.set_figwidth(5.4,forward=True)
    Y=np.sqrt(Yz**2+Yx**2)


    im=ax1.imshow(X,cmap=plt.cm.Blues,extent=[-2,2,-1,1] )

    cuhr=ax1.contour(z,r,np.sqrt(Y**2+X ),np.array([1]),linewidths=2, colors="red"  )
    ax1.clabel(cuhr,inline=True,fmt="UHR",fontsize=12)

    cO=ax1.contour(z,r,X ,np.array([1]),linewidths=2, colors="green"  )
    ax1.clabel(cO,inline=True,fmt="O-cutoff",fontsize=12)


    cR=ax1.contour(z,r,np.sqrt(0.25*Y**2+X )+0.5*Y,np.array([1]),linewidths=2, colors="black"  )
    ax1.clabel(cR,inline=True,fmt='R',fontsize=12)

    cL=ax1.contour(z,r,np.sqrt(0.25*Y**2+X )-0.5*Y,np.array([1]),linewidths=2, colors="black"  )
    ax1.clabel(cL,inline=True,fmt='L-cutoff',fontsize=12)

    ax1.streamplot(z,r,Yz,Yx,minlength=0.9, density=0.4,color="blue")    


    ax1.set_xlabel("z || B_0",fontsize=12)
    ax1.set_ylabel("x",fontsize=12)
    ax1.set_aspect('equal')
    if len(soln)>0:
        plt.plot(soln[:,0],soln[:,1])

    return X
