import numpy as np
import matplotlib.pyplot as plt
from plasma import Plasma



a = 0.38
coils =np.array([-0.63,0.63],np.newaxis) 


p = Plasma(a,coils)

plotX = np.zeros((500,500))
plotY = np.zeros((500,500))
plotYx = np.zeros((500,500))
plotYz = np.zeros((500,500))
plotFlux = np.zeros((500,500))

i  = 0

for z in np.linspace(-0.63,0.63,500):
    j=0
    for x  in np.linspace(-0.38,0.38,500):
    
        plotX[i,j] = p.X(np.array([z,x]))
        plotY[i,j] = p.Yabs(np.array([z,x]))
        plotYx[i,j] = p.Yx(np.array([z,x]))
        plotYz[i,j] = p.Yz(np.array([z,x]))
        plotFlux[i,j] = p.flux(np.array([z,x]))
        j+=1
    i+=1


    #np.where(plotY)


wce=np.where(np.abs(plotY-1)<1e-3)
wce2=np.where(np.abs(plotY-0.5)<1e-3)
wce3=np.where(np.abs(plotY-1./3.)<1e-3)
wce4=np.where(np.abs(plotY-1./4.)<1e-3)
wce5=np.where(np.abs(plotY-1./5.)<1e-3)

wce6=np.where(np.abs(plotY-1./6.)<1e-4)
wce7=np.where(np.abs(plotY-1./7.)<1e-4)
wce8=np.where(np.abs(plotY-1./8.)<1e-4)
wce9=np.where(np.abs(plotY-1./9.)<1e-4)
wce10=np.where(np.abs(plotY-1./10.)<1e-3)
wce11=np.where(np.abs(plotY-1./11.)<1e-4)
wce12=np.where(np.abs(plotY-1./12.)<1e-4)
wce13=np.where(np.abs(plotY-1./13.)<1e-4)
wce14=np.where(np.abs(plotY-1./14.)<1e-3)
wce15=np.where(np.abs(plotY-1./15.)<1e-3)


R = np.where(np.abs(-plotX-plotY+1)<1e-2)
uhr = np.where(np.abs(-plotX-plotY**2+1)<1e-2)

L = np.where(np.abs(-plotX+plotY+1)<1e-2)
O  = np.where(np.abs(-plotX+1)<1e-2)

x = np.linspace(-0.3,0.3,500)
z = np.linspace(-0.6,0.6,500)

    
plt.scatter(z[wce[0]],x[wce[1]],c="green")
plt.scatter(z[R[0]],x[R[1]])
plt.scatter(z[L[0]],x[L[1]],c="g")

plt.scatter(z[O[0]],x[O[1]], c = "r")

plt.scatter(z[wce2[0]],x[wce2[1]], c = "black")
plt.scatter(z[uhr[0]],x[uhr[1]], c = "yellow")
plt.scatter(z[wce3[0]],x[wce3[1]], c = "yellow")


plt.scatter(z[wce4[0]],x[wce4[1]], c = "yellow")
plt.scatter(z[wce5[0]],x[wce5[1]], c = "red")

plt.scatter(z[wce6[0]],x[wce6[1]], c = "yellow")
plt.scatter(z[wce7[0]],x[wce7[1]], c = "yellow")
plt.scatter(z[wce8[0]],x[wce8[1]], c = "yellow")
plt.scatter(z[wce9[0]],x[wce9[1]], c = "yellow")
plt.scatter(z[wce10[0]],x[wce10[1]], c = "yellow")
plt.scatter(z[wce11[0]],x[wce11[1]], c = "yellow")
plt.scatter(z[wce12[0]],x[wce12[1]], c = "yellow")
plt.scatter(z[wce13[0]],x[wce13[1]], c = "yellow")
plt.scatter(z[wce14[0]],x[wce14[1]], c = "yellow")
plt.scatter(z[wce15[0]],x[wce15[1]], c = "yellow")



