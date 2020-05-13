import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys
import os
sys.path.insert(0,'/home/rumiantcev/RiP_general_new')

from conf import p,mf
from dispRel import A,B,C, disp_rel
from coldDispRel import cold_disp_rel,coldNx_Xm,coldDispRel_pol,S,coldNx_Om
import refractive_index as ri
import os
sys.path.insert(0,os.path.abspath("engine/hot"))
import params

import engine as e
import dispRel_EBW as ebw
path = os.path.abspath("results/slab/harm1")
files=os.listdir(os.path.abspath(path))


soln=np.load(path+"/"+files[0]).item()["soln"]
soln_init_cond=np.load(path+"/"+files[0]).item()["init_cond"]
params.Y = float(files[0][files[0].index("_")+1:-4])
NN = []
acc=soln_init_cond[2:4]
dt = 10./2e4
for i in range(0,soln[:,0].size)[0:2]:
    r = soln[i,0:2]
    N = soln[i,2:4]

    X=p.X(r)
    Y=mf.Yabs(r)
    Nz,Nx=N

    acc=acc+e.f(soln[i,0:2].reshape(2,1),soln[i,2:4].reshape(2,1))*dt*200
    #print acc[1]
    NN.append(acc)
    nz,nx =ri.nz(r.reshape(2,1),N.reshape(2,1)),ri.nx(r.reshape(2,1),N.reshape(2,1))
    if np.abs(nx)<1:
        print soln[i,1],disp_rel(X,Y,p.gamma(10),nz,nx),cold_disp_rel(X,Y,nz,nx,"Xm")
    else:
        print soln[i,1],ebw.disp_rel(X,Y,p.gamma(10),nz,nx),cold_disp_rel(X,Y,nz,nx,"Xm")
NN=np.array(NN)

