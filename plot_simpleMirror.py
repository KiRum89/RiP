import matplotlib.pyplot as plt
from scipy.special import zeta
import numpy as np
import my_rc
import os
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import sys
import os
sys.path.insert(0,os.path.abspath(""))
from conf import mf,p

#path = os.path.abspath("results/mirror/simple")
#files=os.listdir(os.path.abspath(path))
##the file names have Y at the end. We willl use to sort the files
"""
Y=[]#stranger be aware, Y has nothing to do with mf.
for f in files:
    if f.find("nz")==-1:
        ind=f.index("_")
        Y.append(float(f[ind+1:-4]))
###
"""
e=1
###find the minimum vgz
#soln=np.load(path+"/"+f).item() #soln is an array
#soln=soln["soln"]
#vgz1 = np.diff(soln[:,0])[-1]#last group velocity
"""
for f in np.array(files)[np.argsort(Y)][::e]:
    soln=np.load(path+"/"+f).item() #soln is an array
    soln=soln["soln"]
    vgz2 = np.diff(soln[:,0])[-1]#last group velocity
    if vgz2<vgz1:
	name = f
###
"""

fig=plt.figure()
fig.set_size_inches(5.4,3.33, forward=True)

print Y[::e],np.min(Y[::e])
inds = np.argsort(Y)
Y = np.array(Y)[inds]
i=0
colormap = cm.viridis
normalize = mcolors.Normalize(vmin=np.min(Y[::e]), vmax=np.max(Y[::e]))

arr =[]
for f in np.array(files)[inds][::e]:


    
    color = colormap(normalize(Y[::e][i]))
    

    soln=np.load(path+"/"+f).item() #soln is an array
    soln=soln["soln"]
    arr.append(soln) 
    """
    if np.abs(soln[-1,1])<10:
        plt.plot(soln[:,0],soln[:,1],lw=2,c=color)

    i+=1 
    """
"""
arr = {"soln":soln,"nz0":-0.03,"ncycle":200}
#np.save(os.path.abspath(path)+"/soln_nz1_-0.03",arr)
"""
