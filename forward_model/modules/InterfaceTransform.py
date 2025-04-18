"""
This file contains all the functions about:
(1) wave vector transforms in a layered medium.
(2) modification of amplidudes based on transmission/reflection coefficients (Fresnel's equations).
"""

import numpy as np

# function 1: process of the single fundamental order wave vector from free space to substrate medium.
# input: theta in JCM thetaphi.
# output: kxvalue, kzvalue under JCM 3D coordinate system.
# output2: diffracted angle under our own definition in radian. there is a negative sign compared with JCM thetaphi definition!
def kvectrans1(theta, nin, nout, knormout):
    thetain = np.radians(-theta)
    sinthetaout=nin*np.sin(thetain)/nout
    angradout=np.arcsin(sinthetaout)
    kxout=knormout*np.sin(angradout)
    kzout=-knormout*np.cos(angradout)
    return kxout, kzout, angradout

# function 2: process of the single fundamental order wave vector from substrate medium to free space.
# input: kxin, kzin under JCM 3D coordinate system.
# output: kxout, kzout under JCM 3D coordinate system.
def kvectrans2(kxin, kzin, nin, nout, knormout):
    thetain=np.atan(-kxin/kzin) # curcial point.
    sinthetaout=nin*np.sin(thetain)/nout
    if abs(sinthetaout)>1: # this case will only triger when nin is larger than nout.
        print("Special Case!!")
        angradout=100
        kxout=0
        kzout=0
    else:
        angradout=np.arcsin(sinthetaout)
        kxout=knormout*np.sin(angradout)
        kzout=-knormout*np.cos(angradout)    
    return kxout, kzout, angradout

# function 3: bulk transforms of all diffracted wave vectors from substrate medium to free space.
# input: kxin, kzin under JCM 3D coordinate system.
# output: kxout, kzout under JCM 3D coordinate system.
def kvectrans3(kxin, kzin, nin, nout, knormout):
    thetain=np.atan(-kxin/kzin) # curcial point.
    sinthetaout=nin*np.sin(thetain)/nout
    if np.max(abs(sinthetaout))<=1:
        indexlist= [i for i, tempvalue in enumerate(sinthetaout) if abs(tempvalue)>1]
        angradout=[]
        kxout=[]
        kzout=[]
        for i in range(len(kxin)):
            if i in indexlist:
                kxout.append(0)
                kzout.append(0)
                angradout.append(100)
            else:
                kxvalue=kxin[i]
                kzvalue=kzin[i]
                thetainvalue=np.atan(-kxvalue/kzvalue)
                sinthetaoutvalue=nin*np.sin(thetainvalue)/nout
                angradoutvalue=np.arcsin(sinthetaoutvalue)
                kxoutvalue=knormout*np.sin(angradoutvalue)
                kzoutvalue=-knormout*np.cos(angradoutvalue)
                kxout.append(kxoutvalue)
                kzout.append(kzoutvalue)
                angradout.append(angradoutvalue)
        kxout=np.array(kxout)
        kzout=np.array(kzout)
        angradout=np.array(angradout)
    else:
        print("Special Case!!")
        angradout=np.arcsin(sinthetaout)
        kxout=knormout*np.sin(angradout)
        kzout=-knormout*np.cos(angradout)
    return kxout, kzout, angradout

# function 4.
# input: kxin, kzin (arrays in the substrate), and nin, nout.
# output: refcoeff and transcoeff arrays (to adjust the amplitudes of the total electric fields in the free space).
# Remark: only valid for P-polarization which is the realistic case.
def FresnelFun(kxin, kzin, nin, nout):
    thetain=np.atan(-kxin/kzin) # curcial point.
    sinthetaout=nin*np.sin(thetain)/nout
    if np.max(abs(sinthetaout))<=1:
        indexlist= [i for i, tempvalue in enumerate(sinthetaout) if abs(tempvalue)>1]
        FresnelarrayR=[]
        FresnelarrayT=[]
        for i in range(len(kxin)):
            if i in indexlist:
                FresnelarrayR.append(0)
                FresnelarrayT.append(0)
            else:
                kxvalue=kxin[i]
                kzvalue=kzin[i]
                thetainvalue=np.atan(-kxvalue/kzvalue) # useful angle 1.
                sinthetaoutvalue=nin*np.sin(thetainvalue)/nout
                angradoutvalue=np.arcsin(sinthetaoutvalue) # useful angle 2.

                FresnelvalueR=(nout*np.cos(thetainvalue)-nin*np.cos(angradoutvalue))/(nout*np.cos(thetainvalue)+nin*np.cos(angradoutvalue))
                FresnelarrayR.append(FresnelvalueR)
                FresnelvalueT=2*nin*np.cos(thetainvalue)/(nout*np.cos(thetainvalue)+nin*np.cos(angradoutvalue))
                FresnelarrayT.append(FresnelvalueT)

        FresnelarrayR=np.array(FresnelarrayR)
        FresnelarrayT=np.array(FresnelarrayT)
    else:
        print("Special case!!")
        FresnelarrayR = np.zeros_like(kxin)
        FresnelarrayT = np.ones_like(kxin)

    return FresnelarrayR, FresnelarrayT



