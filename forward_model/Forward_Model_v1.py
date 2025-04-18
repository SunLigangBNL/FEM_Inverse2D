import sys
import os
#import pathlib
import numpy as np
#import itertools
#import time

jcm_root = "/home/sun2024/JCM_2024"
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
#jcmwave.info()

from forward_model.modules.FEMProcessing import IntensityFEM, GetMatmeta
from forward_model.modules.keys import keys as default_keys

# class ConfigParameter:
#     def __init__(self, pitch, line, height, relpermittivity, refindexi, refindext, lambda0, k0, lambda1, k1, dimqxbranch):
#         self.pitch=pitch
#         self.line=line
#         self.height=height
#         self.relpermittivity=relpermittivity
#         self.refindexi=refindexi
#         self.refindext=refindext
#         self.lambda0=lambda0
#         self.k0=k0
#         self.lambda1=lambda1
#         self.k1=k1
#         self.dimqxbranch=dimqxbranch


#     def __repr__(self):
#         params = ", ".join(f"{key}={value}" for key, value in self.__dict__.items())
#         return f"ConfigParameter({params})"


class ForwardModel:
    def __init__(self, config, Multiplicity=1, NThreads=1):
        self.config=config
        # Start daemon at local machine
        jcmwave.daemon.shutdown()
        jcmwave.daemon.add_workstation(
            Hostname = "localhost",
            NThreads = NThreads,
            Multiplicity = Multiplicity)

    #intref, intrefuncertainty=fm.CompRef(thetaarray,config,diroutput,dimqxbranch)
    def CompRef(self, thetaarray, config, diroutput):
        decayswitch=0
        Qxmat, Qzmat, Intensitymat = IntensityFEM(thetaarray, config, diroutput, decayswitch)
        matmeta=GetMatmeta(thetaarray, config, Qxmat, Qzmat, Intensitymat)

        index1=config.dimqxbranch+1
        index2=2*config.dimqxbranch+1
        Qzdomain=matmeta[:,index1:index2,1]
        intmean=matmeta[:,index1:index2,2]
        intuncertainty=0.01*intmean
        print("Done!!")
        return Qzdomain, intmean, intuncertainty

    def ModelEvaluate(self, keys, working_dir=None):

        thetaarray=np.linspace(-20,30,6)
        dirinput="/home/sun2024/Resource/JCM/Inverse_2D/simulation2a/forward_model/jcm"

        calc_keys = default_keys.copy()
        calc_keys.update(keys)
        
        phi=0
        job_ids = []

        for theta in thetaarray:
            if theta%10==0:
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

        

        #intmodel=np.random.rand(len(thetaarray), 5)
        #print("pseudo results generated.")
        return intmodel




