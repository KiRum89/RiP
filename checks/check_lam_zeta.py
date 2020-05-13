from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


from hotDiTen import HotDiTen
from plasma import Plasma
from local_N import Local_N
from derivForHamEq_EBW2 import *

p = Plasma(0.5, np.array([-1,1]))


#X = p.X()


def check_hdt(soln):

    #ind=np.where(np.isnan(soln[:,3]))[0][0]
    #print ind
    #ind = ind - 1
    ind = soln[:-1,3].size-1
    r = [soln[ind-1,0],soln[ind-1,1]]
    Nz = soln[ind-1,2]
    Nx = soln[ind-1,3]

    globalN = np.array([Nz,Nx])
    locN = Local_N(p.Yvec(r), p.gradY(r),globalN)
    hdt = HotDiTen(p.X(r),p.Yabs(r),p.gamma(r),locN.Nz(),locN.Nx(),20)
    print "lam",hdt.lam()
    print "zeta",hdt.zeta(10)

    print "Z", hdt.Z(10)
    print "Zp", hdt.Zp(10)
    print "Zpp", hdt.Zpp(10)
    print "Zpp_naive", hdt.Zpp_naive(10)


    print "X", p.X(r)
    print "dX_dr", p.dX_dr(r)
   
    print "Yabs", p.Yabs(r)
    print "gradY", p.gradY(r)
    print "dY_dr", p.dY_dr(r)
    

    print "lNx",locN.Nx()
    print "lNz",locN.Nz()
    print "dNz_dN",locN.dNz_dN()
    print "dNx_dN",locN.dNx_dN()
    print "dNx_dr",locN.dNx_dr()
    print "dNz_dr",locN.dNz_dr()



    print "Kxx",hdt.Kxx()
    print "Kyy",hdt.Kyy()
    print "Kzz",hdt.Kzz()
    print "Cxy",hdt.Cxy()
    print "Cyz",hdt.Cyz()


    X = p.X(r)
    Y = p.Yabs(r)
    gamma = p.gamma(r)
    X=p.X(r)

    Y=p.Yabs(r)
    gamma = p.gamma(1)
    
    localN = Local_N(p.Yvec(r),p.gradY(r),globalN)
    
    localNz=localN.Nz()
    localNx=localN.Nx()
    dNz_dN=localN.dNz_dN()
    dNx_dN=localN.dNx_dN()
    dNz_dr=localN.dNz_dr()
    dNx_dr=localN.dNx_dr()
    dX_dr = p.dX_dr(r)
    dY_dr = p.dY_dr(r)
    dgamma_dr = p.dgamma_dr(1)
    
    print "vg",-dD_dN( X,Y,gamma,localNz,localNx,dNz_dN,dNx_dN  ) / dD_dw( X,Y,gamma,localNz,localNx  )


    print "dN_dt",dD_dr(X,Y,gamma,localNz,localNx,dNz_dr,dNx_dr,dX_dr,dY_dr,dgamma_dr) / dD_dw( X,Y,gamma,localNz,localNx  ) 




    
