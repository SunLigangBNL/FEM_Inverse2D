"""
This file contains all functions related to the direct FT method and Born approximation.
"""

import numpy as np
from .utils import coordinate_transform, selectsortfun
#from .source_config import SourceConfig
from .geometry_config import GeometryConfig, Scatterer
#from .material_config import MaterialConfig
from .cdsaxs_config import CDSAXSConfig

def sincfun(x):
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 0:
        return 1.0 if x == 0 else np.sin(x) / x
    else: 
        result = np.empty_like(x)
        zero_mask = (x == 0)
        result[zero_mask] = 1.0
        result[~zero_mask] = np.sin(x[~zero_mask]) / x[~zero_mask]
        return result

def delta_approx(x, mu, sigma):
    result=np.exp(-(x-mu)**2 / (2*sigma**2))
    return result

def polygon_area(x, y):
    y_next = np.roll(y, -1)
    x_next = np.roll(x, -1)
    cross_terms = x * y_next - y * x_next
    return 0.5*np.abs(cross_terms.sum())


"""
Compute Form Factor analytically.
# Input: qxvalue, qzvalue, xarray, yarray.
# Output: FT of a 2D indicator function on a general polygon, as an array.
# Remark: 3 cases are considered here. (1) Qx=0. (2) Qx is not zero but Qz=0. (3) Qx and Qz are not 0.
"""
def FormFacPolygon(qxvalue, qzarray, xarray, yarray):
    N=len(xarray) # N-sided polygon.
    intformarray=np.zeros(len(qzarray), dtype=complex)

    if qxvalue==0:
        print("Singularity case for Qx=0 and Qz=0. Shoelace algorithm applied.")
        polyarea=polygon_area(xarray,yarray)
        intformarray=np.ones(len(qzarray))*polyarea
        print("polygon area:",polyarea)
        #print("current intformarray:",intformarray)
    elif np.any(qzarray==0):
        print("Singularity case: Qx nonzero but Qz=0 occurs.")
    
        # split qzarray into nonzero piece and zero piece.
        zero_idx=np.where(qzarray==0)[0]
        nonzero_mask=(qzarray!=0)
        qzarray1=qzarray[nonzero_mask] # qzarray1 does NOT contain zero.
        # now process the nonzero piece.
        intformarray1=np.zeros(len(qzarray1))
        for i in range(N):
            x1, y1=xarray[i], yarray[i]
            x2, y2=xarray[(i+1)%N], yarray[(i+1)%N]
            if x1==x2:
                #print("Special case: integration along vertical boundary.")
                fac1=(-2)*np.exp((-1j)*qxvalue*x1)/(qxvalue*qzarray1)
                fac2=np.exp((-1j)*qzarray1*y2)-np.exp((-1j)*qzarray1*y1)
            else:
                #print("Integration along general boundary.")
                afac=(y2-y1)/(x2-x1)
                bfac=y1-afac*x1
                fac1=(-1)*(1/qzarray1+2*afac/qxvalue)*np.exp((-1j)*qzarray1*bfac)/(qxvalue+qzarray1*afac)
                fac2=np.exp((-1j)*(qxvalue+qzarray1*afac)*x2)-np.exp((-1j)*(qxvalue+qzarray1*afac)*x1)
            intgamma=fac1*fac2
            intformarray1=intformarray1+intgamma
        intformarray[nonzero_mask]=intformarray1 # update 1.
        # now process the zero piece (one element).
        intformvalue2=0
        if zero_idx.size: 
            for i in range(N):
                x1, y1=xarray[i], yarray[i]
                x2, y2=xarray[(i+1)%N], yarray[(i+1)%N]
                if x1==x2:
                    #print("Special case: integration along vertical boundary.")
                    intgamma=(y2-y1)*2j/qxvalue*np.exp(-1j*qxvalue*x1)
                else:
                    #print("Integration along general boundary.")
                    afac=(y2-y1)/(x2-x1)
                    bfac=y1-afac*x1
                    fac1=afac*(np.exp(-1j*qxvalue*x2)*(1/qxvalue**2+x2*1j/qxvalue)-np.exp(-1j*qxvalue*x1)*(1/qxvalue**2+x1*1j/qxvalue))
                    fac2=(-2*afac/qxvalue**2+1j*bfac/qxvalue)*(np.exp(-1j*qxvalue*x2)-np.exp(-1j*qxvalue*x1))
                    intgamma=fac1+fac2
                intformvalue2=intformvalue2+intgamma
            intformarray[zero_idx[0]]=intformvalue2 # update 2.
            
    else:
        #print("Bulk case: no singularity occurred.")
        for i in range(N):
            x1, y1=xarray[i], yarray[i]
            x2, y2=xarray[(i+1)%N], yarray[(i+1)%N]
            if x1==x2:
                #print("Special case: integration along vertical boundary.")
                fac1=(-2)*np.exp((-1j)*qxvalue*x1)/(qxvalue*qzarray)
                fac2=np.exp((-1j)*qzarray*y2)-np.exp((-1j)*qzarray*y1)
            else:
                #print("Integration along general boundary.")
                afac=(y2-y1)/(x2-x1)
                bfac=y1-afac*x1
                fac1=(-1)*(1/qzarray+2*afac/qxvalue)*np.exp((-1j)*qzarray*bfac)/(qxvalue+qzarray*afac) #******
                fac2=np.exp((-1j)*(qxvalue+qzarray*afac)*x2)-np.exp((-1j)*(qxvalue+qzarray*afac)*x1)
            intgamma=fac1*fac2
            intformarray=intformarray+intgamma  
    return intformarray

"""
Main BA code:
(1) intensity based on differential cross section of Born approximation with general grating profile.
(2) layered medium model.
input: qxvalue, qzarray (SI unit).
input: phiradarray, psiradarray. Both are under the sample frame. Both are exactly the variables in Masa's report.
output: I(qx,qz)/I0 array.
"""
def IntBornPolygon(config, qxvaluein, qzarrayin, thetaarray, psiradarray):
    
    pitch=config.pitch
    E0=config.source.E0
    k0=config.source.k0

    deltaqx=np.pi/pitch # per m.
    sigma=0.033*deltaqx # 3*sigma=0.1*deltaqx.
    mu1=round(qxvaluein/deltaqx)*deltaqx
    mu2=0
    
    # Angle-dependent factor based on Masa's report:
    thetaradarray=np.radians(thetaarray) # this is the EXACT phi angle in Masa's note. Check again on 07/20.
    # psi=2theta: Masa's Fig 2 is EXACTLY same with our own angle definition (negative sign to JCM [theta, phi]).
    anglefac=np.cos(psiradarray)*np.cos(psiradarray)/(np.cos(thetaradarray)*np.cos(psiradarray-thetaradarray))
    c0=anglefac*k0**2/(4*pitch**2)

    inttotal=0    
    scatterlayer=config.geometry.scatterlayer # indices of layers containing scatterers. e.g., [1, 5].

    #for m1 in range(len(scatterlayer)):
    for m1 in range(len(config.geometry.scattercoordmat)):
        # contribution from the scatterers.
        intpart1=0
        epsrb=config.material.layerepsrarray[scatterlayer[m1]]
        for m2 in range(len(config.geometry.scattercoordmat[m1])):
            coeff1=config.material.scatterepsrmat[m1][m2]-epsrb # Attention here!
            xarray=config.geometry.scattercoordmat[m1][m2][:,0] # sample frame.
            zarray=config.geometry.scattercoordmat[m1][m2][:,1] # sample frame.
            intformarray=coeff1*FormFacPolygon(qxvaluein, qzarrayin, xarray, zarray)
            intpart1=intpart1+intformarray
        # contribution from the remaining part in the current scatterer layer.
        if epsrb==1:
            intpart2=0
        else:
            xarray=np.array([-0.5*config.pitch, 0.5*config.pitch, 0.5*config.pitch, -0.5*config.pitch])
            zarray=np.array([config.geometry.layerzarray[scatterlayer[m1]-1],config.geometry.layerzarray[scatterlayer[m1]-1],
                             config.geometry.layerzarray[scatterlayer[m1]],config.geometry.layerzarray[scatterlayer[m1]]])
            coeff2=epsrb-1
            intpart2=coeff2*FormFacPolygon(qxvaluein, qzarrayin, xarray, zarray)
        inttotal=inttotal+intpart1+intpart2

    # integral based on the other layers (to exclude: the top air layer, the bottom air layer, the scatter layers).
    intpart3=0
    indexarray1=np.array(config.geometry.layerindexarray) # all layer index.
    indexarray2=np.array([indexarray1[0],indexarray1[-1]])
    indexarray3=np.array(scatterlayer) # index of layer containing scatterers.
    indexarray4=np.setdiff1d(indexarray1, np.union1d(indexarray2,indexarray3))
    
    for n1 in indexarray4:
        coeff3=config.material.layerepsrarray[n1]-1
        xarray=np.array([-0.5*config.pitch, 0.5*config.pitch, 0.5*config.pitch, -0.5*config.pitch]) # sample frame.
        zarray=np.array([config.geometry.layerzarray[n1-1],config.geometry.layerzarray[n1-1],
                         config.geometry.layerzarray[n1],config.geometry.layerzarray[n1]]) # sample frame.
        intformarray=coeff3*FormFacPolygon(qxvaluein, qzarrayin, xarray, zarray)
        intpart3=intpart3+intformarray

    inttotal=inttotal+intpart3
    intensity=c0*np.abs(inttotal)**2
    return intensity

# Complete the full loop of BA. This is a counterpart of IntensityFEM.
# Input: different keys, config, matmeta from FEM.
# Output: Intensitymat, which will be processed again by GetMeta. No!
def IntensityBA(config, keys, matmeta):
    
    Intensitymatnew=np.zeros((len(config.source.thetaarray),len(config.qxpeakarray)))

    # generate confignew based on keys.
    coordmat, centermat, coordmat2=grating_profile(pitch=keys['pitch'], cd=keys['cd'], h=keys['h'], 
                                                   swa=keys['swa'], rtop=keys['r_top'], rbot=keys['r_bot'], nrsamp=5)

    scattercoordmatnew=[]
    vertexmat1 = coordmat2*1e-9
    vertexmatlist=[vertexmat1] # JCM frame.
    scattercoordmatnew.append(coordinate_transform(vertexmatlist)) # sample frame.
    geometrynew=GeometryConfig(scattercoordmat=scattercoordmatnew, layerzarray=config.geometry.layerzarray,
                               scatterlayer=config.geometry.scatterlayer)

    confignew=CDSAXSConfig(pitch=config.pitch, dimqxbranch=config.dimqxbranch, qxpeakarray=config.qxpeakarray, source=config.source, 
                            geometry=geometrynew, material=config.material)

    
    qxindexarray=config.centerindex+config.qxpeakarray
    qxrefarray=np.linspace(-config.dimqxbranch*config.deltaqx,config.dimqxbranch*config.deltaqx,2*config.dimqxbranch+1)
    for i in qxindexarray:
        qxvalue = qxrefarray[i]
        matslice = matmeta[:, i, :]

        Qzarray = matslice[:, 1]
        #intarray = matslice[:, 2]
        thetaarray = matslice[:, 3]
        psiarray = matslice[:, 4]
        flagarray = matslice[:, 5]

        IntBAarray2a=np.zeros(len(Qzarray))

        Qzarray2=Qzarray[flagarray==1]
        #intarray2=intarray[flagarray==1]
        thetaarray2=thetaarray[flagarray==1]
        psiarray2=psiarray[flagarray==1]

        # Now, IntBAarray is corresponding to Qzarray2 format.
        IntBAarray2=IntBornPolygon(confignew, qxvalue, Qzarray2, thetaarray2, psiarray2) # Note new geometry info is used here.

        IntBAarray2a[flagarray==1]=IntBAarray2 # now it has the same length and structure as Qzarray.
        Intensitymatnew[:,i-config.centerindex-1]=IntBAarray2a

    return Intensitymatnew

def grating_profile(pitch, cd, h, swa, rtop, rbot, nrsamp):
    
    x1=-pitch/2
    
    shift1=h/(2*np.tan(np.radians(swa)))
    shift2=h/(2*np.tan(np.radians(swa)))
    x2=-cd/2-shift1
    x3=-cd/2+shift1
    x4=cd/2-shift2
    x5=cd/2+shift2

    x6=pitch/2

    # original coordinates of the trapezoid.
    coordmat=np.array([[x1, 0],
                       [x2, 0],
                       [x3, h],
                       [x4, h],
                       [x5, 0],
                       [x6, 0]])
    
    vecmat=coordmat[1:]-coordmat[:-1] # 5 by 2. vectors of each boundary.
    norms=np.linalg.norm(vecmat, axis=1).reshape(-1, 1)
    uvecmat=vecmat/norms # 5 by 2. unit vectors of each boundary.

    # coordinates of the curve surface.
    coordmat2=np.zeros((nrsamp*4+2,2))
    coordmat2[0,:]=coordmat[0,:]
    coordmat2[-1,:]=coordmat[-1,:]

    centermat=np.zeros((4,2))
    
    for i in range(1,5): # always generate 4 circles.
        #print("i=",i)
        # compute arc angle.
        dotproduct=np.dot(uvecmat[i-1,:], uvecmat[i,:])
        arcangle=np.arccos(dotproduct) # always positive.
        #print("arc angle in degree (always less than 180) : ", np.degrees(arcangle))

        if i == 1:
            r=rbot #****************
            theta1=-np.pi/2 #****************
            theta2=-np.pi/2+arcangle #****************
        elif i==2:
            r=rtop #****************
            theta1=np.pi/2+arcangle #****************
            theta2=np.pi/2 #****************
        elif i==3:
            r=rtop #****************
            theta1=np.pi/2 #****************
            theta2=np.pi/2-arcangle #****************
        else:
            r=rbot #****************
            theta1=-np.pi/2-arcangle #****************
            theta2=-np.pi/2 #****************  
        
        #print("current radius : ",r)
        anglearray=np.linspace(theta1,theta2,nrsamp)

        # compute center of the incircle or excircle.
        normvec=uvecmat[i,:]-uvecmat[i-1,:] # center of the incircle or excircle is located on this line. Attention: NOT a unit vector here!
        tempangle=np.arccos(np.dot(normvec,uvecmat[i,:])/(np.linalg.norm(normvec))) # corrected version.
        #print("tempangle in degree (always less than 90) : ",np.degrees(tempangle))
        deltax=r/np.tan(tempangle) # x-coordinate of center shifted a bit. always positive.
        #print("deltax : ",deltax)

        if i==1:
            center=np.array([coordmat[i,0]-deltax,0+r]) #*******************
        elif i==2:
            center=np.array([coordmat[i,0]+deltax,h-r]) #*******************
        elif i==3:
            center=np.array([coordmat[i,0]-deltax,h-r]) #*******************
        else:
            center=np.array([coordmat[i,0]+deltax,0+r]) #*******************
        #print("center :",center)
        centermat[i-1,:]=center

        # coordinates of the sampling points on the arc.
        xarray=center[0]+r*np.cos(anglearray)
        yarray=center[1]+r*np.sin(anglearray)

        index1=1+(i-1)*nrsamp
        index2=index1+nrsamp
        coordmat2[index1:index2,0] = xarray
        coordmat2[index1:index2,1] = yarray
    return coordmat, centermat, coordmat2


# Version 3: introduce non-symmetry for both rtop/rbot and SWA.
def grating_profile3(pitch, cd, h, swaleft, swaright, rtopleft, rtopright, rbotleft, rbotright, nrsamp):
    
    x1=-pitch/2
    
    shift1=h/(2*np.tan(np.radians(swaleft)))
    shift2=h/(2*np.tan(np.radians(swaright)))
    x2=-cd/2-shift1
    x3=-cd/2+shift1
    x4=cd/2-shift2
    x5=cd/2+shift2

    x6=pitch/2

    # original coordinates of the trapezoid.
    coordmat=np.array([[x1, 0],
                       [x2, 0],
                       [x3, h],
                       [x4, h],
                       [x5, 0],
                       [x6, 0]])
    
    vecmat=coordmat[1:]-coordmat[:-1] # 5 by 2. vectors of each boundary.
    norms=np.linalg.norm(vecmat, axis=1).reshape(-1, 1)
    uvecmat=vecmat/norms # 5 by 2. unit vectors of each boundary.

    # coordinates of the curve surface.
    coordmat2=np.zeros((nrsamp*4+2,2))
    coordmat2[0,:]=coordmat[0,:]
    coordmat2[-1,:]=coordmat[-1,:]

    centermat=np.zeros((4,2))
    
    for i in range(1,5): # always generate 4 circles.
        # compute arc angle.
        dotproduct=np.dot(uvecmat[i-1,:], uvecmat[i,:])
        arcangle=np.arccos(dotproduct) # always positive.

        if i == 1:
            r=rbotleft #****************
            theta1=-np.pi/2 #****************
            theta2=-np.pi/2+arcangle #****************
        elif i==2:
            r=rtopleft #****************
            theta1=np.pi/2+arcangle #****************
            theta2=np.pi/2 #****************
        elif i==3:
            r=rtopright #****************
            theta1=np.pi/2 #****************
            theta2=np.pi/2-arcangle #****************
        else:
            r=rbotright #****************
            theta1=-np.pi/2-arcangle #****************
            theta2=-np.pi/2 #****************  
        
        #print("current radius : ",r)
        anglearray=np.linspace(theta1,theta2,nrsamp)

        # compute center of the incircle or excircle.
        normvec=uvecmat[i,:]-uvecmat[i-1,:] # center of the incircle or excircle is located on this line. Attention: NOT a unit vector here!
        tempangle=np.arccos(np.dot(normvec,uvecmat[i,:])/(np.linalg.norm(normvec))) # corrected version.
        deltax=r/np.tan(tempangle) # x-coordinate of center shifted a bit. always positive.

        if i==1:
            center=np.array([coordmat[i,0]-deltax,0+r]) #*******************
        elif i==2:
            center=np.array([coordmat[i,0]+deltax,h-r]) #*******************
        elif i==3:
            center=np.array([coordmat[i,0]-deltax,h-r]) #*******************
        else:
            center=np.array([coordmat[i,0]+deltax,0+r]) #*******************
        centermat[i-1,:]=center

        # coordinates of the sampling points on the arc.
        xarray=center[0]+r*np.cos(anglearray)
        yarray=center[1]+r*np.sin(anglearray)

        index1=1+(i-1)*nrsamp
        index2=index1+nrsamp
        coordmat2[index1:index2,0] = xarray
        coordmat2[index1:index2,1] = yarray
    return coordmat, centermat, coordmat2
    
    
# # Far field intensity based on BA: only valid for rectangular geometry.
# # input: qxvalue, qzarray. (unit in per angstrom).
# # output: I(qx,qz) array.
# # Remark: this is the correct version based on Masa's note.
# def IntBornRec(config, qxvaluein, qzarrayin, thetaarray, psiradarray):

#     # unit transformation.
#     qxvalue=qxvaluein*10**10 # now unit is per m.
#     qzarray=qzarrayin*10**10 # now unit is per m.
    
#     pitch=config.pitch
#     line=config.line
#     height=config.height
#     H=config.thickness
#     E0=config.E0
#     k0=config.k0
#     rr=config.sdd
#     epsilonr1=config.epsilonr1
#     epsilonr2=config.epsilonr2

#     deltaqx=np.pi/pitch # per m.
#     sigma=0.033*deltaqx # 3*sigma=0.1*deltaqx.

#     ninc=config.refindexi
#     ntrm=config.refindext

#     # Angle-dependent factor introduced by Masa:
#     thetaradarray=np.radians(thetaarray) # this is the EXACT phi angle in Masa's note.
#     psiradarray=np.array(psiradarray) # this is the EXACT psi angle in Masa's note. LHS kf is cooresponding to a positive psi angle.
#     anglefac=np.cos(psiradarray)*np.cos(psiradarray)/(np.cos(thetaradarray)*np.cos(psiradarray-thetaradarray))

#     c0=anglefac*k0**2/(4*pitch**2)

#     # compute the kernel integral.
#     c1=1-np.real(epsilonr1)
#     c2=1-np.real(epsilonr2)
#     # contribution from the grating strucutre.
#     mu1=round(qxvalue/deltaqx)*deltaqx
#     c3=c1*line*height*sincfun(qxvalue*line/2)*sincfun(qzarray*height/2)*np.exp(1j*qzarray*height/2)*delta_approx(qxvalue,mu1,sigma)
#     # contribution from the substrate in [-H,0].
#     mu2=0
#     c4=c2*2*np.pi*H*sincfun(qzarray*H/2)*np.exp(1j*qzarray*H/2)*delta_approx(qxvalue,mu2,sigma)
#     inttemp=np.abs(c3+c4)
#     intensity=c0*inttemp**2

#     return intensity

