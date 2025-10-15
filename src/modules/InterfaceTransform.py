"""
This file contains the following functions:
(1) wave vector transform from final substrate layer into free space.
(2) modification of amplidudes based on transmission/reflection coefficients (Fresnel's equations).
"""

import numpy as np

"""
function 1: Transform wave vectors from substrate medium to free space.
Inputs: (JCM 3D frame)
    kxarrayin, kzarrayin: wave vectors under JCM 3D coordinate system.
    refindexin: refractive index of input medium.
    refindexout: refractive index of output medium (usually 1).
    knormout: |k| in output medium = 2π/λ * n_out.
Returns:
    kxarrayout, kzarrayout (JCM 3D frame). 
    angradout: diffraction angles in radian. (With a negative sign to the JCM [theta, phi] definition.)
"""
def karraytrans(kxarrayin, kzarrayin, refindexin, refindexout, knormout):
    if (np.real(refindexin) > 1.0) or (np.abs(np.imag(refindexin)) > 0):
        print(f"ERROR: refective index of the substratge = {refindexin} is >1 or has an imaginary part. Aborted.")
        sys.exit(1)
        
    thetain = np.arctan2(kxarrayin, -kzarrayin) # from JCM angle to our defined angle frame.
    sinthetaout = refindexin * np.sin(thetain) / refindexout
    angradarrayout=np.arcsin(sinthetaout)
    kxarrayout=knormout*np.sin(angradarrayout)
    kzarrayout=-knormout*np.cos(angradarrayout)
    return kxarrayout, kzarrayout, angradarrayout

"""
function 2: transmission coefficients.
input: kxin, kzin (arrays in the substrate), and nin, nout.
output: refcoeff and transcoeff arrays (to adjust the amplitudes of the total electric fields in the free space).
Remark: only valid for P-polarization, which is the realistic case.
"""
def FresnelFun(kxarrayin, kzarrayin, refindexin, refindexout):
    if (np.real(refindexin) > 1.0) or (np.abs(np.imag(refindexin)) > 0):
        print(f"ERROR: refective index of the substratge = {refindexin} is >1 or has an imaginary part. Aborted.")
        sys.exit(1)
    
    thetain = np.arctan2(kxarrayin, -kzarrayin) # from JCM angle to our defined angle frame.      
    sinthetaout = refindexin*np.sin(thetain)/refindexout           
    angradout = np.arcsin(sinthetaout)   # principal branch
    cosin  = np.cos(thetain)
    cosout = np.cos(angradout)                 
    temp = refindexout*cosin + refindexin*cosout  
    # P‑polarized Fresnel formulas
    FresnelarrayR = (refindexout*cosin-refindexin*cosout)/temp
    FresnelarrayT = 2*refindexin*cosin/temp
    return FresnelarrayR, FresnelarrayT

