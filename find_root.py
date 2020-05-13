import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from dispRel import disp_rel, A,B,C
#from hotDiTen import HotDiTen
from plasma import SlabPlasma, Plasma, PlasmaParabMF


from coldDispRel import cold_disp_rel, coldDispRel_pol, coldNx_Xm, coldNx_Om

coils = np.array([-0.63,0.63])

#p = Plasma(0.38,coils)
#p = PlasmaParabMF(1)
p = SlabPlasma(1)

"""
print "hot, not a soluion",disp_rel(p.X([0,1.01]),p.Yabs([0,1.01]),p.gamma(1),0.6,0.8)
print "cold, not a soluion",coldDispRel_pol(p.X([0,1.01]),p.Yabs([0,1.01]),0.6,0.8)


print cold_disp_rel(p.X([0,1.01]),p.Yabs([0,1.01]),0.6,0.8,"Xm")## disp relation with branches explicitly 

coldNxXm=np.sqrt(coldNx_Xm(p.X([0,1.01]),p.Yabs([0,1.01]),0.6))
print "cold Nx**2 Xm", coldNxXm
print "hot disp rel for cold Nx Xm", disp_rel(p.X([0,1.01]),p.Yabs([0,1.01]),p.gamma(1),0.6,coldNxXm) 
print "cold disp rel, no branches",coldDispRel_pol(p.X([0,1.01]),p.Yabs([0,1.01]),0.6,coldNxXm) 

coldNxOm=np.sqrt(coldNx_Om(p.X([0,1.01]),p.Yabs([0,1.01]),0.6))
print "cold Nx**2 Om", coldNxOm
print "hot disp rel for cold Nx Om", disp_rel(p.X([0,1.01]),p.Yabs([0,1.01]),p.gamma(1),0.6,coldNxOm) 





print "cold, AAH, Xm" ,cold_disp_rel(p.X([0,1.01]),p.Yabs([0,1.01]),0.6,coldNxXm,"Xm")


print "cold, AAH, Om" ,cold_disp_rel(p.X([0,1.01]),p.Yabs([0,1.01]),0.6,coldNxOm,"Om")
"""


def get_ABC(X,Y,gamma,Nx,Nz):
    from hotDiTen import HotDiTen
    from dispRel import A,B,C
    hdt = HotDiTen(X,Y,gamma,Nx,Nz,10)
    dielTen = [hdt.Kxx(),hdt.Kyy(),hdt.Kzz(),hdt.Cxy(),hdt.Cxz(),hdt.Cyz()]
    return  A(Nz, dielTen),B(Nz, dielTen),C(Nz, dielTen)
    


coldNxXm=np.sqrt(coldNx_Xm(p.X([0,1.139]),p.Yabs([0,1.139]),0.6))
A,B,C=get_ABC(p.X([0,1.139]),p.Yabs([0,1.139]),p.gamma(1),0.6,coldNxXm)

print coldNxXm
print B**2 - 4*A*C 




i = 0
for x in np.linspace(0,5,1000):
    Nx_Xm=coldNx_Xm(p.X([0,x]),p.Yabs([0,x]),0.01)

    Nx_Om=coldNx_Om(p.X([0,x]),p.Yabs([0,x]),0.01)


    print Nx_Om,Nx_Xm
    #plt.scatter(p.X([0,x]),Nx_Xm)
    #plt.scatter(p.X([0,x]),Nx_Om)
    plt.scatter(x,Nx_Xm)
    plt.scatter(x,Nx_Om)
 
    i +=1


x = -0.18
Nx_Xm=coldNx_Xm(p.X([0,x]),p.Yabs([0,x]),0.1)**2

Nx_Om=coldNx_Om(p.X([0,x]),p.Yabs([0,x]),0.1)**2



#plt.ylim([0,1000])

