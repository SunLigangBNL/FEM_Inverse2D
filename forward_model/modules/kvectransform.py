#! /usr/bin/env python3

""" Functions to transform wave vectors between different medium. """

import numpy as np


# function 1: given theta (JCM value), compute the transmitted wavevector under fundamental order.
# input: theta in JCM thetaphi.
# output: kxvalue, kzvalue under JCM coordinate.
# output2: diffracted angle under our own definition in radian. there is a negative sign compared with jcm thetaphi definition!
# This is a better version based on kvectrans1.
def kvectrans1a(theta, nin, nout, knormout):
    thetain = np.radians(-theta)
    sinthetaout=nin*np.sin(thetain)/nout
    angradout=np.arcsin(sinthetaout)
    kxout=knormout*np.sin(angradout)
    kzout=-knormout*np.cos(angradout)
    return kxout, kzout, angradout

# function 2.
# input: kxin, kzin being scalars.
# output: kxout, kzout, angradout being scalars.
def kvectrans2(kxin, kzin, nin, nout, knormout):
    thetain=np.atan(-kxin/kzin) # curcial point.
    sinthetaout=nin*np.sin(thetain)/nout
    if abs(sinthetaout)>1: # this case will never be triggered when nin is smaller than nout.
        angradout=100
        kxout=0
        kzout=0
    else:
        angradout=np.arcsin(sinthetaout)
        kxout=knormout*np.sin(angradout)
        kzout=-knormout*np.cos(angradout)    
    return kxout, kzout, angradout

# function 3.
# input: kxin, kzin being arrays.
# output: kxout, kzout, angradout being arrays.
# Remark: this function need to be modified when nin is larger than nout.
def kvectrans2a(kxin, kzin, nin, nout, knormout):
    thetain=np.atan(-kxin/kzin) # curcial point.
    sinthetaout=nin*np.sin(thetain)/nout
    if np.max(abs(sinthetaout))>1:
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
        angradout=np.arcsin(sinthetaout)
        kxout=knormout*np.sin(angradout)
        kzout=-knormout*np.cos(angradout)
    return kxout, kzout, angradout