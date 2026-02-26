"""
This file contains all functions related to the direct FT method and Born approximation.
"""

import numpy as np
from .utils import coordinate_transform, selectsortfun, sortfunBA
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
        #print("Singularity case for Qx=0 and Qz=0. Shoelace algorithm applied.")
        polyarea=polygon_area(xarray,yarray)
        intformarray=np.ones(len(qzarray))*polyarea
    elif np.any(qzarray==0):
        #print("Singularity case: Qx nonzero but Qz=0 occurs.")
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
def IntBornPolygon(config, qxvaluein, qzarrayin):
    
    pitch=config.pitch
    k0=config.source.k0
    
    # Angle-dependent factor based on Masa's report:
    # thetaradarray=np.radians(thetaarray) # this is the EXACT phi angle in Masa's note. Check again on 07/20.
    # psi=2theta: Masa's Fig 2 is EXACTLY same with our own angle definition (negative sign to JCM [theta, phi]).
    # anglefac=np.cos(psiradarray)*np.cos(psiradarray)/(np.cos(thetaradarray)*np.cos(psiradarray-thetaradarray))

    # # alternative method: compute anglefac by using Qx, Qz directly.
    # Qlength=2*k0*np.sin(np.abs(psiradarray/2))
    # tempfac1=(Qlength/(2*k0))**2
    # anglefac2=(1-2*tempfac1)**2/((qxvaluein/Qlength)**2-tempfac1)

    # alternative method: compute anglefac by using Qx, Qz directly.
    Qlength=np.sqrt(qxvaluein**2+qzarrayin**2)
    tempfac1=(Qlength/(2*k0))**2
    anglefac2=(1-2*tempfac1)**2/((qxvaluein/Qlength)**2-tempfac1)

    # temporarily turn off the angle-dependent factor.
    # anglefac2=1.0
    
    c0=anglefac2*k0**2/(4*pitch**2)

    inttotal=0    
    scatterlayer=config.geometry.scatterlayer # indices of layers containing scatterers. e.g., [1, 5].

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
# Input: different keys, config, externally computed Qzmat.
# Output: Intensitymat.
# Remark: the Qzgridmat is computed independently without using FEM matmeta.
# Remark2: geometry flipped for NIST case.

def IntensityBA(keys, config, Qzgridmat):
    thetaarray=config.source.thetaarray
    qxrefarray=np.linspace(-config.dimqxbranch*config.deltaqx,config.dimqxbranch*config.deltaqx,2*config.dimqxbranch+1)
    Intensitymatnew=np.zeros((len(thetaarray),len(qxrefarray)))
    
    # rewrite profile coordinates for different rounding models.
    if config.geometrymode==0:
        scattercoordmatnew=config.geometry.scattercoordmat # default polygon model with JCM coordinates. already flipped if needed.
    elif config.geometrymode==1:
        print("Rounding model detected: Type 1.")
        scattercoordmatnew=[]
        orgmat, centermat, surfmat=SurfaceCoordinates1(keys)
        if config.flipswitch==1:
            surfmat2=surfmat.copy()
            surfmat2[:, 1]=-surfmat[:,1]
            surfmatsorted=surfmat2[1:-1]
        else:
            surfmat2=surfmat[::-1,:]
            surfmat3=surfmat2[1:-1]
            surfmatsorted=np.vstack([surfmat3[-1:], surfmat3[:-1]])
        vertexmat1=surfmatsorted*1e-9    
        vertexmatlist=[vertexmat1] # JCM frame.
        scattercoordmatnew.append(coordinate_transform(vertexmatlist)) # sample frame.
            
    elif config.geometrymode==2:
        print("Rounding model detected: Type 2.")
        scattercoordmatnew=[]
        orgmat, centermat, surfmat=SurfaceCoordinates2(keys)
        if config.flipswitch==1:
            surfmat2=surfmat.copy()
            surfmat2[:, 1]=-surfmat[:,1]
            surfmatsorted=surfmat2[1:-1]
        else:
            surfmat2=surfmat[::-1,:]
            surfmat3=surfmat2[1:-1]
            surfmatsorted=np.vstack([surfmat3[-1:], surfmat3[:-1]])
        vertexmat1=surfmatsorted*1e-9    
        vertexmatlist=[vertexmat1] # JCM frame.
        scattercoordmatnew.append(coordinate_transform(vertexmatlist)) # sample frame.
        
    elif config.geometrymode==3:
        print("Rounding model detected: Type 3.")
        scattercoordmatnew=[]
        orgmat, centermat, surfmat=SurfaceCoordinates3(keys)
        if config.flipswitch==1:
            surfmat2=surfmat.copy()
            surfmat2[:, 1]=-surfmat[:,1]
            surfmatsorted=surfmat2[1:-1]
        else:
            surfmat2=surfmat[::-1,:]
            surfmat3=surfmat2[1:-1]
            surfmatsorted=np.vstack([surfmat3[-1:], surfmat3[:-1]])
        vertexmat1=surfmatsorted*1e-9    
        vertexmatlist=[vertexmat1] # JCM frame.
        scattercoordmatnew.append(coordinate_transform(vertexmatlist)) # sample frame.
    
    elif config.geometrymode==4:
        print("Rounding model detected: Type 4.")
        scattercoordmatnew=[]
        orgmat, centermat, surfmat=SurfaceCoordinates3(keys)
        
        if config.flipswitch==1:
            surfmat2=surfmat.copy()
            surfmat2[:, 1]=-surfmat[:,1]
            surfmatsorted=surfmat2[1:-1]
            
            #index1=keys['Narcsamp']
            #index2=index1+keys['Ntr']

            index1 = int(np.array(keys['Narcsamp']).item())
            index2 = index1 + int(np.array(keys['Ntr']).item())

            N = surfmatsorted.shape[0]
            j1 = (N - 1) - index1   # mirror of index1
            j2 = (N - 1) - index2   # mirror of index2
            
            mat1 = np.vstack([surfmatsorted[:index1+1], surfmatsorted[j1:]])
            mat2 = np.vstack([surfmatsorted[index1:index2+1], surfmatsorted[j2:j1+1]])
            mat3 = surfmatsorted[index2:j2+1].copy()
            
            vertexmat1=mat1*1e-9
            vertexmat2=mat2*1e-9
            vertexmat3=mat3*1e-9
            
        else:
            surfmat2=surfmat[::-1,:]
            surfmat3=surfmat2[1:-1]
            surfmatsorted=np.vstack([surfmat3[-1:], surfmat3[:-1]])
            
            # index1=keys['Narcsamp']+1
            # index2=index1+keys['Ntr']
            
            index1 = int(np.array(keys['Narcsamp']).item())+1
            index2 = index1 + int(np.array(keys['Ntr']).item())
            
            N = surfmatsorted.shape[0]
            y1 = surfmatsorted[index1, 1]
            y2 = surfmatsorted[index2, 1]
            idx1_pair = [i for i in np.where(np.isclose(surfmatsorted[:, 1], y1))[0] if i != index1][0]
            idx2_pair = [i for i in np.where(np.isclose(surfmatsorted[:, 1], y2))[0] if i != index2][0]
            if idx1_pair < index1:
                idx1_pair = [i for i in np.where(np.isclose(surfmatsorted[:, 1], y1))[0] if i != index1][-1]
            if idx2_pair < index2:
                idx2_pair = [i for i in np.where(np.isclose(surfmatsorted[:, 1], y2))[0] if i != index2][-1]
            mat4 = np.vstack([surfmatsorted[:index1+1], surfmatsorted[idx1_pair:]])
            mat5 = np.vstack([surfmatsorted[index1:index2+1], surfmatsorted[idx2_pair:idx1_pair+1]])
            mat5 = np.vstack([mat5[-1:], mat5[:-1]])
            mat6 = surfmatsorted[index2:idx2_pair+1].copy()
            mat6 = np.vstack([mat6[-1:], mat6[:-1]])

            vertexmat1=mat4*1e-9
            vertexmat2=mat5*1e-9
            vertexmat3=mat6*1e-9
        
        vertexmatlist=[vertexmat1, vertexmat2, vertexmat3] # JCM frame.***************************************** new model here!
        scattercoordmatnew.append(coordinate_transform(vertexmatlist)) # sample frame.
    else:
        raise ValueError("Parameter config.geometrymode is NOT set correctly.")
    
    geometrynew=GeometryConfig(scattercoordmat=scattercoordmatnew, layerzarray=config.geometry.layerzarray,
                               scatterlayer=config.geometry.scatterlayer)
    confignew=CDSAXSConfig(pitch=config.pitch, dimqxbranch=config.dimqxbranch, qxpeakarray=config.qxpeakarray,
                           obindexmask=config.obindexmask, flipswitch=config.flipswitch, geometrymode=config.geometrymode,
                           source=config.source, geometry=geometrynew, material=config.material)

    for i in range(len(qxrefarray)): # for all standard dimensions.
        Qxvalue=qxrefarray[i]
        Qzarray=Qzgridmat[:,i]
        mask = Qzarray != 0 # only pass non-zero elements to the IntBornPolygon function.
        Qzarray2=Qzarray[mask]
        IntBAarray2=IntBornPolygon(confignew, Qxvalue, Qzarray2)
        Intensitymatnew[mask, i]=IntBAarray2

    return Intensitymatnew


def GetQzgrid(config):
    thetaarray=config.source.thetaarray
    qxrefarray=np.linspace(-config.dimqxbranch*config.deltaqx,config.dimqxbranch*config.deltaqx,2*config.dimqxbranch+1)
    Qzgridmat=np.zeros((len(thetaarray),len(qxrefarray)))
    
    tempfac1=config.source.lambda0/config.pitch
    
    QxmatBA=[]
    QzmatBA=[]
    for i in range(len(thetaarray)):
        theta=thetaarray[i]
        thetarad=np.radians(theta)
        diffangrad0=-thetarad # under sample frame.
        kivecx4=config.source.k0*np.sin(diffangrad0) # under sample frame.
        kivecz4=config.source.k0*np.cos(diffangrad0) # under sample frame.
    
        # compute diffradarray based on diffraction equation.
        tempfac2=np.sin(diffangrad0)
        m1=(1-tempfac2)/tempfac1
        m2=(-1-tempfac2)/tempfac1
        M1=np.ceil(min(m1,m2))
        M2=np.floor(max(m1,m2))
    
        indexarray=np.arange(M1,M2+1)
        sindiffradarray=tempfac2+indexarray*tempfac1
        cosdiffradarray=np.abs(np.sqrt(1-sindiffradarray**2))
    
        kxarray4=config.source.k0*sindiffradarray
        kzarray4=config.source.k0*cosdiffradarray
        Qxarray=kxarray4-kivecx4
        Qzarray=kzarray4-kivecz4
    
        indexlist1=np.argsort(Qxarray)
        Qxarray2=Qxarray[indexlist1]
        Qzarray2=Qzarray[indexlist1]
    
        QxmatBA.append(Qxarray2)
        QzmatBA.append(Qzarray2)
        
    QxmatBA=np.array(QxmatBA, dtype=object)
    QzmatBA=np.array(QzmatBA, dtype=object)

    matmetaBA=np.zeros((len(thetaarray),2*config.dimqxbranch+1,4)) # per angle, per Qx location, store (Qx, Qz, theta, flag) values.
    for i in range(len(thetaarray)):
        qxarray=np.array(QxmatBA[i])
        for j,qxvalue in enumerate(qxrefarray):
            index1 = np.argmin(np.abs(qxarray-qxvalue))
            if np.abs(qxarray[index1]-qxvalue)<0.5*config.deltaqx: # the 0.5 here is good for both periodic grating and finite grating.
                matmetaBA[i,j,0]=QxmatBA[i][index1]
                matmetaBA[i,j,1]=QzmatBA[i][index1]
                matmetaBA[i,j,2]=thetaarray[i]
                if QzmatBA[i][index1]==0 and thetaarray[i]==0:
                    matmetaBA[i,j,3]=0
                else:
                    matmetaBA[i,j,3]=1
    
    for i in range(len(qxrefarray)):
        matsliceBA=matmetaBA[:, i, :]
        QzarrayBA=sortfunBA(matsliceBA)
        Qzgridmat[:,i]=QzarrayBA
    
    return Qzgridmat

# Compute coordinates of the surface for Model1 (symmetric single trapezoid with rounded corners).
# Input: keys.
# Output: trapezoid coordinates, centers of all circle of curvature, curved surface coordinates.
def SurfaceCoordinates1(keys):
    pitch=keys['pitch']
    cd=keys['cd']
    h=keys['h']
    swa=keys['swa']
    rtop=keys['rtop']
    rbot=keys['rbot']
    nrsamp=keys['nrsamp']
    
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
        # compute arc angle.
        dotproduct=np.dot(uvecmat[i-1,:], uvecmat[i,:])
        arcangle=np.arccos(dotproduct) # always positive.

        if i == 1:
            r=rbot
            theta1=-np.pi/2 
            theta2=-np.pi/2+arcangle 
        elif i==2:
            r=rtop 
            theta1=np.pi/2+arcangle 
            theta2=np.pi/2 
        elif i==3:
            r=rtop 
            theta1=np.pi/2 
            theta2=np.pi/2-arcangle 
        else:
            r=rbot 
            theta1=-np.pi/2-arcangle 
            theta2=-np.pi/2 

        anglearray=np.linspace(theta1,theta2,nrsamp)

        # compute center of the incircle or excircle.
        normvec=uvecmat[i,:]-uvecmat[i-1,:] # center of the incircle or excircle is located on this line. Attention: NOT a unit vector here!
        tempangle=np.arccos(np.dot(normvec,uvecmat[i,:])/(np.linalg.norm(normvec))) # corrected version.
        deltax=r/np.tan(tempangle) # x-coordinate of center shifted a bit. always positive.

        if i==1:
            center=np.array([coordmat[i,0]-deltax,0+r])
        elif i==2:
            center=np.array([coordmat[i,0]+deltax,h-r]) 
        elif i==3:
            center=np.array([coordmat[i,0]-deltax,h-r])
        else:
            center=np.array([coordmat[i,0]+deltax,0+r])
        centermat[i-1,:]=center

        # coordinates of the sampling points on the arc.
        xarray=center[0]+r*np.cos(anglearray)
        yarray=center[1]+r*np.sin(anglearray)

        index1=1+(i-1)*nrsamp
        index2=index1+nrsamp
        coordmat2[index1:index2,0] = xarray
        coordmat2[index1:index2,1] = yarray
    return coordmat, centermat, coordmat2

# Compute coordinates of the surface for Model2 (Asymmetric single trapezoid with rounded corners).
# Input: keys.
# Output: trapezoid coordinates, centers of all circle of curvature, curved surface coordinates.
def SurfaceCoordinates2(keys):
    pitch=keys['pitch']
    cd=keys['cd']
    h=keys['h']
    swaleft=keys['swaleft']
    swaright=keys['swaright']
    rtopleft=keys['rtopleft']
    rtopright=keys['rtopright']
    rbotleft=keys['rbotleft']
    rbotright=keys['rbotright']
    nrsamp=keys['nrsamp']
    
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
    
    for i in range(1,5):
        # compute arc angle.
        dotproduct=np.dot(uvecmat[i-1,:], uvecmat[i,:])
        arcangle=np.arccos(dotproduct) # always positive.

        if i == 1:
            r=rbotleft 
            theta1=-np.pi/2 
            theta2=-np.pi/2+arcangle 
        elif i==2:
            r=rtopleft 
            theta1=np.pi/2+arcangle 
            theta2=np.pi/2 
        elif i==3:
            r=rtopright
            theta1=np.pi/2
            theta2=np.pi/2-arcangle
        else:
            r=rbotright
            theta1=-np.pi/2-arcangle
            theta2=-np.pi/2
        
        anglearray=np.linspace(theta1,theta2,nrsamp)

        # compute center of the incircle or excircle.
        normvec=uvecmat[i,:]-uvecmat[i-1,:] # center of the incircle or excircle is located on this line. Attention: NOT a unit vector here!
        tempangle=np.arccos(np.dot(normvec,uvecmat[i,:])/(np.linalg.norm(normvec))) # corrected version.
        deltax=r/np.tan(tempangle) # x-coordinate of center shifted a bit. always positive.

        if i==1:
            center=np.array([coordmat[i,0]-deltax,0+r])
        elif i==2:
            center=np.array([coordmat[i,0]+deltax,h-r])
        elif i==3:
            center=np.array([coordmat[i,0]-deltax,h-r])
        else:
            center=np.array([coordmat[i,0]+deltax,0+r])
        centermat[i-1,:]=center

        # coordinates of the sampling points on the arc.
        xarray=center[0]+r*np.cos(anglearray)
        yarray=center[1]+r*np.sin(anglearray)

        index1=1+(i-1)*nrsamp
        index2=index1+nrsamp
        coordmat2[index1:index2,0] = xarray
        coordmat2[index1:index2,1] = yarray
    return coordmat, centermat, coordmat2

# Compute coordinates of the surface for Model3 (Asymmetric stacked trapezoids with rounded corners).
# Input: keys.
# Output: trapezoid coordinates, centers of all circle of curvature, curved surface coordinates.
# example keys:

def SurfaceCoordinates3(keys):
    pitch=keys['pitch']
    Narcsamp=int(round(keys['Narcsamp']))
    Ntr=int(round(keys['Ntr']))
    h=keys['htot']
    heightb=keys['hbot']
    heightt=keys['htop']
    swabl=keys['swabl']
    swatl=keys['swatl']
    swabr=keys['swabr']
    swatr=keys['swatr']
    rbl=keys['rbl']
    rtl=keys['rtl']
    rbr=keys['rbr']
    rtr=keys['rtr']
    xl0=keys['xl0']
    xr0=keys['xr0']

    Ntemp=Ntr+4 # number of total vertices without rounding. 
    Ntotal=2*(Ntr+1)+4*Narcsamp+2 # number of total vertices with rounding.
    
    coordmat=np.zeros((2*Ntemp,2))
    centermat=np.zeros((4,2))
    coordmat2=np.zeros((Ntotal,2))
    
    harray=np.zeros(Ntr+4)
    harray[2]=heightb
    harray[Ntr+3]=h
    harray[2:Ntr+3]=np.linspace(heightb, h-heightt, Ntr+1)
    
    xarrayl=np.zeros(Ntr+4)
    xarrayr=np.zeros(Ntr+4)
    
    xarrayl[0]=-pitch/2
    xarrayl[2]=xl0
    shiftbl=heightb/(np.tan(np.radians(swabl)))
    xarrayl[1]=xarrayl[2]-shiftbl

    if abs(xarrayl[1])>=pitch/2:    
        raise ValueError("shiftbl exceeds the unit cell!!")
    
    xarrayr[0]=pitch/2
    xarrayr[2]=xr0
    shiftbr=heightb/(np.tan(np.radians(swabr)))
    xarrayr[1]=xarrayr[2]+shiftbr

    if abs(xarrayr[1])>=pitch/2:    
        raise ValueError("shiftbr exceeds the unit cell!!")
    
    for i in range(1, Ntr+1):
        xarrayl[i+2]=keys[f'xl{i}']
        xarrayr[i+2]=keys[f'xr{i}']
    
    shifttl=heightt/(np.tan(np.radians(swatl)))
    xarrayl[-1]=xarrayl[-2]+shifttl
    
    shifttr=heightt/(np.tan(np.radians(swatr)))
    xarrayr[-1]=xarrayr[-2]-shifttr
    
    coordmat[:,0]=np.concatenate((xarrayl, xarrayr[::-1]))
    coordmat[:,1]=np.concatenate((harray, harray[::-1]))

    indexarray=np.array([0,Ntemp-2,Ntemp-1,2*Ntemp-3])
    indexarray2=np.array([1, 1+Narcsamp*1+Ntr+1, 1+Narcsamp*2+Ntr+1, Ntotal-Narcsamp-1])
    
    vecmat=coordmat[1:]-coordmat[:-1] # vectors of each boundary.
    norms=np.linalg.norm(vecmat, axis=1).reshape(-1, 1)
    uvecmat=vecmat/norms # unit vectors of each boundary.
    
    coordmat2=np.zeros((Ntotal,2))
    coordmat2[0,:]=coordmat[0,:]
    coordmat2[-1,:]=coordmat[-1,:]
    coordmat2[indexarray2[0]+Narcsamp:indexarray2[1],:]=coordmat[indexarray[0]+1+1:indexarray[0]+1+1+Ntr+1,:]
    coordmat2[indexarray2[2]+Narcsamp:indexarray2[3],:]=coordmat[indexarray[2]+1+1:indexarray[2]+1+1+Ntr+1,:]

    # check the 4 radiuses.
    thresholdbl=heightb/np.sin(np.radians(swabl))*np.tan(np.radians((180-swabl)/2))
    thresholdbr=heightb/np.sin(np.radians(swabr))*np.tan(np.radians((180-swabr)/2))
    thresholdtl=heightt/np.sin(np.radians(swatl))*np.tan(np.radians((180-swatl)/2))
    thresholdtr=heightt/np.sin(np.radians(swatr))*np.tan(np.radians((180-swatr)/2))

    if rbl>thresholdbl:
        #print("Parameter rbl is too large !!")
        raise ValueError("Parameter rbl is too large.")
    
    if rbr>thresholdbr:
        #print("Parameter rbr is too large !!")
        raise ValueError("Parameter rbr is too large.")

    if rtl>thresholdtl:
        #print("Parameter rtl is too large !!")
        raise ValueError("Parameter rtl is too large.")

    if rtr>thresholdtr:
        #print("Parameter rtr is too large !!")
        raise ValueError("Parameter rtr is too large.")
    
    for i in range(4): # always generate 4 circles.
        index=indexarray[i]
        dotproduct=np.dot(uvecmat[index,:], uvecmat[index+1,:])
        arcangle=np.arccos(dotproduct) # arc angle.
        if i==0:
            r=rbl
            theta1=-np.pi/2
            theta2=-np.pi/2+arcangle
        elif i==1:
            r=rtl 
            theta1=np.pi/2+arcangle
            theta2=np.pi/2 
        elif i==2:
            r=rtr 
            theta1=np.pi/2 
            theta2=np.pi/2-arcangle
        else:
            r=rbr 
            theta1=-np.pi/2-arcangle 
            theta2=-np.pi/2
        # compute center of the incircle or excircle.
        anglearray=np.linspace(theta1,theta2,Narcsamp)
        normvec=uvecmat[index+1,:]-uvecmat[index,:] # center of the incircle/excircle is located on this line. NOT a unit vector here!
        tempangle=np.arccos(np.dot(normvec,uvecmat[index+1,:])/(np.linalg.norm(normvec))) # corrected version.
        deltax=r/np.tan(tempangle) # x-coordinate of center shifted.
        
        if i==0:
            center=np.array([coordmat[index+1,0]-deltax,0+r])
        elif i==1:
            center=np.array([coordmat[index+1,0]+deltax,h-r])
        elif i==2:
            center=np.array([coordmat[index+1,0]-deltax,h-r]) 
        else:
            center=np.array([coordmat[index+1,0]+deltax,0+r]) 
        centermat[i,:]=center
    
        # coordinates of the sampling points on the arc.
        xarray=center[0]+r*np.cos(anglearray)
        yarray=center[1]+r*np.sin(anglearray)

        if max(abs(xarray))>=pitch/2:
            #print("Some arc exceeds the unit cell!!")
            raise ValueError("Some arc exceeds the unit cell!!")
        
        index2=indexarray2[i]
        coordmat2[index2:index2+Narcsamp,0]=xarray
        coordmat2[index2:index2+Narcsamp,1]=yarray
    return coordmat, centermat, coordmat2  
