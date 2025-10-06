"""
This file contains the functions transforming the FEM solution in the k space to the Q space.
"""

import os
import numpy as np

from .InterfaceTransform import karraytrans, FresnelFun

# function1: transform FEM solution to intensity in Q space.
# input: thetaarray, JCM Output folder, decayswitch, sample detector distance (unit: m).
# output: Qxmat, Qzmat, Intensitymat. (unit: per angstrom.)
# What is exactly contained in Intensitymat? I_{t,m}/I_i=Tm^2|km,z|/(E0^2*k_{i,z}). See the note at page D6.3.

def IntensityFEM(config, directory):
    Qxmat=[]
    Qzmat=[]
    Intensitymat=[]
    Psiradmat=[]
    thetaarray=config.source.thetaarray
    for i in range(len(thetaarray)):
        theta=thetaarray[i]
        thetarad=np.radians(theta)
        # if i%20==0:
        #     print("current thata : ",theta)    
        foldername = os.path.join(directory, f"theta_{theta:.2f}")
        pathkx=os.path.join(foldername, 'kxarray.txt')
        pathky=os.path.join(foldername, 'kyarray.txt')
        pathkz=os.path.join(foldername, 'kzarray.txt')
        
        pathexspec=os.path.join(foldername, 'exspec.txt')
        patheyspec=os.path.join(foldername, 'eyspec.txt')
        pathezspec=os.path.join(foldername, 'ezspec.txt')

        # Transform coordinates: from JCM 2D to JCM 3D:
        kxarray = np.loadtxt(pathkx,dtype=complex)
        kzarray = np.loadtxt(pathky,dtype=complex) # note the changed name.
        kxarray=np.real(kxarray)
        kzarray=np.real(kzarray)
        
        exspec = np.genfromtxt(pathexspec, dtype=complex)
        eyspec = np.genfromtxt(pathezspec, dtype=complex) # note the changed name.
        ezspec = np.genfromtxt(patheyspec, dtype=complex) # note the changed name.

        # Release diffracted waves from substrate into free space.
        # layerrefindexarray[-1] is always air.
        refindexlastlayer=np.real(config.material.layerrefindexarray[-2]) # Must take the real part to coinside with the JCM setting.
        refindexair=1.0
        kxarray3, kzarray3, diffradarray = karraytrans(kxarrayin=kxarray, kzarrayin=kzarray, refindexin=refindexlastlayer, refindexout=refindexair,
                                                       knormout=config.source.k0) 

        # transform diffradarray to psiarray. (psiarray is the angle respecting to the transmitted kivec!!)
        diffangrad0=-thetarad # diffraction angle of the fundamental order wave kivec (with a negative sign to the JCM [theta, phi] definition).
        kivecx0=config.source.k0*np.sin(diffangrad0) # JCM 3D frame
        kivecz0=-config.source.k0*np.cos(diffangrad0) # JCM 3D frame
        psiradarray=diffradarray-diffangrad0 # psi=2theta. sample frame.
        
        FresnelR, FresnelT=FresnelFun(kxarrayin=kxarray,kzarrayin=kzarray,refindexin=refindexlastlayer, refindexout=refindexair)
        
        # modify the amplitudes of the transmitted plane waves.       
        exspec2=FresnelT*exspec
        eyspec2=FresnelT*eyspec
        ezspec2=FresnelT*ezspec
        
        # compute intensity array based on diffraction efficiency and conservation of energy.
        temparray1=np.real(exspec2*np.conjugate(exspec2)+eyspec2*np.conjugate(eyspec2)+ezspec2*np.conjugate(ezspec2))
        temparray2=np.abs(kzarray3) # here we should use free-space wave vector.
        intarray=temparray1*temparray2/(config.source.E0**2*config.source.k0*np.abs(np.cos(thetarad)))

        # coordinate transformation: from JCM 3D frame to CDSAXS sample frame.
        kxarray4=kxarray3
        kzarray4=-kzarray3
        kivecx4=kivecx0
        kivecz4=-kivecz0
        
        Qxarray=kxarray4-kivecx4
        Qzarray=kzarray4-kivecz4
        
        # sort the datasets based on increasing Qxvalue.
        indexlist1=np.argsort(Qxarray)
        Qxarray2=Qxarray[indexlist1]
        Qzarray2=Qzarray[indexlist1]
        intarray2=intarray[indexlist1]
        psiradarray2=psiradarray[indexlist1]
        
        Qxmat.append(Qxarray2)
        Qzmat.append(Qzarray2)
        Intensitymat.append(intarray2)
        Psiradmat.append(psiradarray2)
    
    Qxmat=np.array(Qxmat, dtype=object)
    Qzmat=np.array(Qzmat, dtype=object)
    Intensitymat=np.array(Intensitymat, dtype=object)
    Psiradmat=np.array(Psiradmat, dtype=object)

    return Qxmat, Qzmat, Intensitymat, Psiradmat

# function2:
# input: Qxmat, Qzmat, Intensitymat (SI unit).
# output: matmeta in regular grid (SI unit).
# Remark: matmeta[:,:,3] contains the original JCM thetaarray. (not the defined sample rotation angle.)
# Remark: matmeta[:,:,4] contains delta diffraction angle 2theta.
def GetMatmeta(config, Qxmat, Qzmat, Intensitymat, Psiradmat):
    thetaarray=config.source.thetaarray
    deltaqx=config.deltaqx
    qxrefarray=np.linspace(-config.dimqxbranch*config.deltaqx,config.dimqxbranch*config.deltaqx,2*config.dimqxbranch+1)
    
    matmeta=np.zeros((len(thetaarray),2*config.dimqxbranch+1,6)) # per angle, per Qx location, store (Qx, Qz, Intensity, theta, psi, flag) values.
    for i in range(len(thetaarray)):
        qxarray=np.array(Qxmat[i])
        for j,qxvalue in enumerate(qxrefarray):
            index1 = np.argmin(np.abs(qxarray-qxvalue))
            if np.abs(qxarray[index1]-qxvalue)<0.5*deltaqx: # the 0.5 here is good for both periodic grating and finite grating.
                matmeta[i,j,0]=Qxmat[i][index1]
                matmeta[i,j,1]=Qzmat[i][index1]
                matmeta[i,j,2]=Intensitymat[i][index1]
                matmeta[i,j,3]=thetaarray[i]
                matmeta[i,j,4]=Psiradmat[i][index1]
                if Qzmat[i][index1]==0 and Intensitymat[i][index1]==0 and thetaarray[i]==0 and Psiradmat[i][index1]==0:
                    matmeta[i,j,5]=0
                else:
                    matmeta[i,j,5]=1
                
    print("matmeta generated.")
    return matmeta
