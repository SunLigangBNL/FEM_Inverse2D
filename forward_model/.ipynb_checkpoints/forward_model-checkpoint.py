#! /usr/bin/env python3

import sys
import os
import pathlib
import numpy as np
import itertools
import time

jcm_root = "/home/sun2024/JCM_2024" # -> set JCMROOT installation directory.
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
#jcmwave.info()

from forward_model.modules.kvectransform import kvectrans1a, kvectrans2, kvectrans2a
from forward_model.modules.keys import keys as default_keys
#from forward_model.modules.load_data import load_data

class ConfigParameter:
    def __init__(self, pitch, line, height, relpermittivity, refindexi, refindext, lambda0, k0, lambda1, k1, dimqxbranch):
        self.pitch=pitch
        self.line=line
        self.height=height
        self.relpermittivity=relpermittivity
        self.refindexi=refindexi
        self.refindext=refindext
        self.lambda0=lambda0
        self.k0=k0
        self.lambda1=lambda1
        self.k1=k1
        self.dimqxbranch=dimqxbranch


    def __repr__(self):
        params = ", ".join(f"{key}={value}" for key, value in self.__dict__.items())
        return f"ConfigParameter({params})"


class ForwardModel:
    def __init__(self, config, NThreads=1, Multiplicity=1):
        #self._data = load_data()
        self.config=config
        # Start daemon at local machine
        jcmwave.daemon.shutdown()
        jcmwave.daemon.add_workstation(
            Hostname = "localhost",
            NThreads = NThreads,
            Multiplicity = Multiplicity)

    # Compute numerical reference. Corresponding to the exp_data() function in JCM benchmark.
    # Input: thetaarray.
    # Output: intensity_mean, intensity_uncertainty
    def CompRef(self, thetaarray, config, diroutput, dimqxtemp):
        Qxmat=[]
        Qzmat=[]
        Intensitymat=[]
        
        ninc=config.refindexi
        ntrm=config.refindext
        k0=config.k0
        k1=config.k1

        for i in range(len(thetaarray)):
            theta=thetaarray[i]
            thetarad=np.radians(theta)
                
            foldername = os.path.join(diroutput, f"theta_{theta:.2f}")
            pathkx=os.path.join(foldername, 'kxarray.txt')
            pathky=os.path.join(foldername, 'kyarray.txt')
            #pathkz=os.path.join(foldername, 'kzarray.txt')
            
            pathexspec=os.path.join(foldername, 'exspec.txt')
            patheyspec=os.path.join(foldername, 'eyspec.txt')
            pathezspec=os.path.join(foldername, 'ezspec.txt')
            
            kxarray = np.loadtxt(pathkx,dtype=complex)
            kzarray = np.loadtxt(pathky,dtype=complex) # note the changed name.
            kxarray=np.real(kxarray)
            kzarray=np.real(kzarray)
            
            exspec = np.genfromtxt(pathexspec, dtype=complex)
            eyspec = np.genfromtxt(pathezspec, dtype=complex) # note the changed name.
            ezspec = np.genfromtxt(patheyspec, dtype=complex) # note the changed name.
        
            intarray=np.real(exspec*np.conjugate(exspec)+eyspec*np.conjugate(eyspec)+ezspec*np.conjugate(ezspec))
        
            # apply the kvectrans functions.
            kivecx2, kivecz2, diffangrad2 = kvectrans1a(theta, ninc, ntrm, k1)
            kivecx3, kivecz3, diffangrad3 = kvectrans2(kivecx2, kivecz2, ntrm, ninc, k0) # note the order of the parameters!!
            kxarray3, kzarray3, anglearray = kvectrans2a(kxarray, kzarray, ntrm, ninc, k0)
            
            # coordinate transformation.
            kqarrayx=-kxarray3
            kqarrayz=kzarray3
            
            kivecqx=-kivecx3
            kivecqz=kivecz3
            
            Qxarray=kqarrayx-kivecqx # unit: per meter.
            Qzarray=kqarrayz-kivecqz # unit: per meter.
            Qxarray2=10**(-10)*Qxarray # transform unit to per angstrom.
            Qzarray2=10**(-10)*Qzarray # transform unit to per angstrom.

            indexlist1=np.argsort(Qxarray2)
            Qxarray3=Qxarray2[indexlist1]
            Qzarray3=Qzarray2[indexlist1]
            intarray3=intarray[indexlist1]
            
            Qxmat.append(Qxarray3)
            Qzmat.append(Qzarray3)
            Intensitymat.append(intarray3)
        
        Qxmat=np.array(Qxmat, dtype=object)
        Qzmat=np.array(Qzmat, dtype=object)
        Intensitymat=np.array(Intensitymat, dtype=object)
        print("Intensity in Q-space computed.")
        
        periodicityx=config.pitch*10**10 # unit: angstrom.
        deltaqx=2*np.pi/periodicityx # unit: per angstrom.
        dimqx=2*dimqxtemp+1 # number of all interested Qx.
        print("number of interested Qx values in total: ", dimqx)
        qxrefarray=np.linspace(-dimqxtemp*deltaqx,dimqxtemp*deltaqx,dimqx)
        
        matmeta=np.zeros((len(thetaarray),dimqx,3)) # per angle, per Qx location, store (Qx, Qz, Intensity) values.
        for i in range(len(thetaarray)):
            qxarray=np.array(Qxmat[i])
            for j,qxvalue in enumerate(qxrefarray):
                index1 = np.argmin(np.abs(qxarray-qxvalue))
                if np.abs(qxarray[index1]-qxvalue)<0.5*deltaqx:
                    matmeta[i,j,0]=Qxmat[i][index1]
                    matmeta[i,j,1]=Qzmat[i][index1]
                    matmeta[i,j,2]=Intensitymat[i][index1]
        
        Qzmat=matmeta[:,:,1]
        int_mean=matmeta[:,:,2]
        print("all expectation computed.")
        int_uncertainty=int_mean*np.random.uniform(0.3,0.6)
        print("all uncertainty generated.")
        return Qzmat, int_mean, int_uncertainty
    
    def ModelEvaluate(self, keys, working_dir=None):
        # for now, we locally introduce config and dirinput and thetaarray.

        thetaarray=np.linspace(-20,30,6)
        
        dirinput="/home/sun2024/Resource/JCM/Inverse_2D/simulation1c/forward_model/jcm2"

        # geometry parameter: unit m.
        pitch=80e-9
        
        # material parameter.
        relpermittivity=0.999435 #**************update
        refindexi=1.0
        refindext=np.sqrt(relpermittivity)
        
        # source definition: unit m.
        lambda0=1e-9 #********************************************** make sure this is correct!!
        k0=2*np.pi/lambda0
        lambda1=lambda0/refindext
        k1=2*np.pi/lambda1
        
        dimqxbranch=2 # number of interested Qx peak per branch.

        ninc=refindexi
        ntrm=refindext

        # Qxmat=[]
        # Qzmat=[]
        # Intensitymat=[]
 
        calc_keys = default_keys.copy()
        calc_keys.update(keys)
        
        phi=0
        job_ids = []

        for theta in thetaarray:
            if theta%6==0:
                print("current thata : ",theta)
            these_calc_keys = calc_keys.copy()
            these_calc_keys["theta"] = theta
            these_calc_keys["phi"] = phi

            if working_dir is not None:
                wd = os.path.join(working_dir,'phi{}_theta{}'.format(phi,theta))
            else:
                wd = None
                
            job_id = jcmwave.solve(os.path.join(dirinput, "project.jcmpt"),
                                   keys=these_calc_keys,
                                   temporary=(working_dir is None),
                                   working_dir=wd)
            job_ids.append(job_id)

        # wait for the calculations
        job_statuses = jcmwave.daemon.status(job_ids)
        print("All status",job_statuses)
        results, logs = jcmwave.daemon.wait(job_ids=job_ids, verbose=False)

        # export the log file.
        log_filename = os.path.join(dirinput, "logs_all.txt")
        with open(log_filename, "w") as log_file:
            sys.stdout = log_file
            for i, log in enumerate(logs):
                print(f"Log {i}:")
                print(log["Log"]["Out"])
                print("*" * 100)
            sys.stdout = sys.__stdout__
        print(f"All logs saved in {log_filename}.")

        intmean=np.random.rand(len(thetaarray), 5)
        print("pseudo results generated.")
        return intmean




