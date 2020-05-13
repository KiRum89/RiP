import sys
import os
sys.path.insert(0,os.path.abspath("plasma"))
sys.path.insert(0,os.path.abspath("magneticField"))
sys.path.insert(0,os.path.abspath(""))
sys.path.insert(0,os.path.abspath("diTen/cold"))
sys.path.insert(0,os.path.abspath("diTen/hot"))
sys.path.insert(0,os.path.abspath("plotting"))
sys.path.insert(0,os.path.abspath("checks"))
sys.path.insert(0,os.path.abspath("engine/hot"))
sys.path.insert(0,os.path.abspath("engine/cold"))
sys.path.insert(0,os.path.abspath("scenarios"))



magnetic_field_model ='slab'


plasma_model = 'slab'

verbose=True #to control some printed message

if plasma_model== "slab":
    p = __import__("plasma_slab")

if plasma_model == "mirror":

    p = __import__("plasma")


if plasma_model == "mirror_bvp":

    p = __import__("plasma_bvp")




if magnetic_field_model == "slab":
    mf = __import__("mf_slab")
    


    
if magnetic_field_model == "mirror":
    
    I1, I2, I = 500,500,320
    #I1, I2, I = 1000,1000,200
    name = str(I1)+"_"+str(I2)+"_"+str(I)

    mf = __import__("mf_mirror")

if magnetic_field_model == "simple":
    mf = __import__("mf_mirrorSimple")

print plasma_model
print magnetic_field_model
