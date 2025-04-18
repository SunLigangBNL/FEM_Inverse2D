"""
This file contains the functions transforming the FEM solution in the k space to the Q space.
"""

import os
import numpy as np

from .InterfaceTransform import kvectrans1, kvectrans2, kvectrans3, FresnelFun

# function1: transform FEM solution to intensity in Q space.
# input: thetaarray, JCM Output folder, decayswitch, sample detector distance (unit: m).
# output: Qxmat, Qzmat, Intensitymat. (unit: per angstrom.)
def IntensityFEM(thetaarray, config, folderoutput, decayswitch):
    Qxmat=[]
    Qzmat=[]
    Intensitymat=[]
    for i in range(len(thetaarray)):
        theta=thetaarray[i]
        thetarad=np.radians(theta)
        if i%20==0:
            print("current thata : ",theta)    
        foldername = os.path.join(folderoutput, f"theta_{theta:.2f}")
        pathkx=os.path.join(foldername, 'kxarray.txt')
        pathky=os.path.join(foldername, 'kyarray.txt')
        #pathkz=os.path.join(foldername, 'kzarray.txt')
        
        pathexspec=os.path.join(foldername, 'exspec.txt')
        patheyspec=os.path.join(foldername, 'eyspec.txt')
        pathezspec=os.path.join(foldername, 'ezspec.txt')

        # Transform coordinates: from JCM 2D to JCM 3D:
        kxarray = np.loadtxt(pathkx,dtype=complex)
        kzarray = np.loadtxt(pathky,dtype=complex)
        kxarray=np.real(kxarray)
        kzarray=np.real(kzarray)
        
        exspec = np.genfromtxt(pathexspec, dtype=complex)
        eyspec = np.genfromtxt(pathezspec, dtype=complex) # note the changed name.
        ezspec = np.genfromtxt(patheyspec, dtype=complex) # note the changed name.
    
    
        # apply interface transform functions.
        ninc=config.refindexi
        ntrm=config.refindext
        kivecx2, kivecz2, diffangrad2 = kvectrans1(theta, ninc, ntrm, config.k1)
        kivecx3, kivecz3, diffangrad3 = kvectrans2(kivecx2, kivecz2, ntrm, ninc, config.k0) # note the order of the parameters!!
        kxarray3, kzarray3, anglearray = kvectrans3(kxarray, kzarray, ntrm, ninc, config.k0) # note that we need to use the modified one!
        FresnelR, FresnelT=FresnelFun(kxarray,kzarray,ntrm, ninc)
    
        # modify the amplitudes of the transmitted plane waves.       
        exspec2=FresnelT*exspec
        eyspec2=FresnelT*eyspec
        ezspec2=FresnelT*ezspec
        
        # compute the far-field intensity.
        intarray=np.real(exspec2*np.conjugate(exspec2)+eyspec2*np.conjugate(eyspec2)+ezspec2*np.conjugate(ezspec2))
        
        # coordinate transformation: from JCM 3D to CDSAXS.
        kqarrayx=-kxarray3
        kqarrayz=kzarray3
        
        kivecqx=-kivecx3
        kivecqz=kivecz3
        
        Qxarray=kqarrayx-kivecqx
        Qzarray=kqarrayz-kivecqz
        
        Qxarray2=10**(-10)*Qxarray
        Qzarray2=10**(-10)*Qzarray
        
        # sort the datasets based on increasing Qxvalue.
        indexlist1=np.argsort(Qxarray2)
        Qxarray3=Qxarray2[indexlist1]
        Qzarray3=Qzarray2[indexlist1]
        
        intarray2=intarray[indexlist1] # now intarray2 has the sorted order.
        anglearray2=anglearray[indexlist1] # now anglearray2 has the sorted order.
        #print("angle array in FEM",anglearray2)
        if decayswitch==1:
            #print("Decay factor 1/r included.")
            distancearray=config.rr/np.cos(anglearray2)
            decayarray=1/distancearray
            intarray3=decayarray*intarray2
        elif decayswitch==2:
            #print("Decay factor 1/r^2 included.")
            distancearray=config.rr/np.cos(anglearray2)
            decayarray=1/distancearray**2
            intarray3=decayarray*intarray2
        else:
            #print("No decay factor included.")
            intarray3=intarray2
        
        Qxmat.append(Qxarray3)
        Qzmat.append(Qzarray3)
        Intensitymat.append(intarray3)
    
    Qxmat=np.array(Qxmat, dtype=object)
    Qzmat=np.array(Qzmat, dtype=object)
    Intensitymat=np.array(Intensitymat, dtype=object)

    return Qxmat, Qzmat, Intensitymat

# function2:
# input: Qxmat, Qzmat, Intensitymat (unit per angstrom).
# output: matmeta in regular grid (unit per angstrom).
def GetMatmeta(thetaarray, config, Qxmat, Qzmat, Intensitymat):

    periodicityx=config.pitch*10**10 # unit: angstrom.
    deltaqx=2*np.pi/periodicityx # unit: per angstrom.
    dimqx=2*config.dimqxbranch+1 # number of all interested Qx.
    print("number of interested Qx values in total: ", dimqx)
    qxrefarray=np.linspace(-config.dimqxbranch*deltaqx,config.dimqxbranch*deltaqx,dimqx)
    
    matmeta=np.zeros((len(thetaarray),dimqx,4)) # per angle, per Qx location, store (Qx, Qz, Intensity, theta) values.
    for i in range(len(thetaarray)):
        qxarray=np.array(Qxmat[i])
        for j,qxvalue in enumerate(qxrefarray):
            index1 = np.argmin(np.abs(qxarray-qxvalue))
            if np.abs(qxarray[index1]-qxvalue)<0.5*deltaqx:
                matmeta[i,j,0]=Qxmat[i][index1]
                matmeta[i,j,1]=Qzmat[i][index1]
                matmeta[i,j,2]=Intensitymat[i][index1]
                matmeta[i,j,3]=thetaarray[i]
    print("matmeta generated.")
    return matmeta

# function3: transform results from dameon (instead of directory) to intensity in Q space.
# This function will only be used inside ModelEvaluate.
# input: thetaarray, dameon results, decayswitch.
# output: Qxmat, Qzmat, Intensitymat. (unit: per angstrom.)
def IntensityFEM2(thetaarray, config, results, decayswitch):
    Qxmat=[]
    Qzmat=[]
    Intensitymat=[]
    for i in range(len(thetaarray)):
        theta=thetaarray[i]
        
        result=results[i]
        farresult=result[3] # corresponding to the "fourier_modes_minus_y.jcm" in the project.jcmp file.
        karray=farresult["K"]
        earray=farresult["ElectricFieldStrength"][0]


        kxarray=np.real(karray[:,0])
        kzarray=np.real(karray[:,1]) # note the changed name.
    
        exspec=earray[:,0]
        ezspec=earray[:,1] # note the changed name.
        eyspec=earray[:,2] # note the changed name.
    
        # apply interface transform functions.
        ninc=config.refindexi
        ntrm=config.refindext
        kivecx2, kivecz2, diffangrad2 = kvectrans1(theta, ninc, ntrm, config.k1)
        kivecx3, kivecz3, diffangrad3 = kvectrans2(kivecx2, kivecz2, ntrm, ninc, config.k0) # note the order of the parameters!!
        kxarray3, kzarray3, anglearray = kvectrans3(kxarray, kzarray, ntrm, ninc, config.k0) # note that we need to use the modified one!
        FresnelR, FresnelT=FresnelFun(kxarray,kzarray,ntrm, ninc)
    
        # modify the amplitudes of the transmitted plane waves.      
        exspec2=FresnelT*exspec
        eyspec2=FresnelT*eyspec
        ezspec2=FresnelT*ezspec
        
        # compute the far-field intensity.
        intarray=np.real(exspec2*np.conjugate(exspec2)+eyspec2*np.conjugate(eyspec2)+ezspec2*np.conjugate(ezspec2))
        
        # coordinate transformation: from JCM 3D to CDSAXS.
        kqarrayx=-kxarray3
        kqarrayz=kzarray3
        
        kivecqx=-kivecx3
        kivecqz=kivecz3
        
        Qxarray=kqarrayx-kivecqx
        Qzarray=kqarrayz-kivecqz
        
        Qxarray2=10**(-10)*Qxarray
        Qzarray2=10**(-10)*Qzarray
        
        # sort the datasets based on increasing Qxvalue.
        indexlist1=np.argsort(Qxarray2)
        Qxarray3=Qxarray2[indexlist1]
        Qzarray3=Qzarray2[indexlist1]
        
        intarray2=intarray[indexlist1] # now intarray2 has the sorted order.
        anglearray2=anglearray[indexlist1] # now anglearray2 has the sorted order.
        if decayswitch==1:
            #print("Decay factor 1/r included.")
            distancearray=config.rr/np.cos(anglearray2)
            decayarray=1/distancearray
            intarray3=decayarray*intarray2
        elif decayswitch==2:
            #print("Decay factor 1/r^2 included.")
            distancearray=config.rr/np.cos(anglearray2)
            decayarray=1/distancearray**2
            intarray3=decayarray*intarray2
        else:
            #print("No decay factor included.")
            intarray3=intarray2
        
        Qxmat.append(Qxarray3)
        Qzmat.append(Qzarray3)
        Intensitymat.append(intarray3)
    
    Qxmat=np.array(Qxmat, dtype=object)
    Qzmat=np.array(Qzmat, dtype=object)
    Intensitymat=np.array(Intensitymat, dtype=object)

    return Qxmat, Qzmat, Intensitymat
