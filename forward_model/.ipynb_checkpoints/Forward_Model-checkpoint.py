import sys
import os
import numpy as np

jcm_root = "/home/sun2024/JCM_2024"
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
#jcmwave.info()

from forward_model.modules.FEMProcessing import IntensityFEM, IntensityFEM2, GetMatmeta
from forward_model.modules.keys import keys as default_keys

class ForwardModel:
    def __init__(self, config, Multiplicity=1, NThreads=1):
        self.config=config
        # Start daemon at local machine
        jcmwave.daemon.shutdown()
        jcmwave.daemon.add_workstation(
            Hostname = "localhost",
            NThreads = NThreads,
            Multiplicity = Multiplicity)

    # Compute numerical reference.
    def CompRef(self, thetaarray, config, diroutput):
        decayswitch=0
        Qxmat, Qzmat, Intensitymat = IntensityFEM(thetaarray, config, diroutput, decayswitch)
        matmeta=GetMatmeta(thetaarray, config, Qxmat, Qzmat, Intensitymat)

        # select from the first peak.
        index1=config.dimqxbranch+1
        index2=2*config.dimqxbranch+1
        Qzdomain=matmeta[:,index1:index2,1]
        intmean=matmeta[:,index1:index2,2]
        intuncertainty=intmean*np.random.uniform(0.3,0.6)
        print("Numerical reference generated.")
        return Qzdomain, intmean, intuncertainty

    # Main forward model.
    # Input: single key (fixed problem setup).
    # Process: rotation scan of multiple simulations.
    # Output: corresponding intensity map.
    def ModelEvaluate(self, keys, config, working_dir=None):

        # preparation.
        thetaarray=np.linspace(-20,30,6)
        dirinput="/home/sun2024/Resource/JCM/Inverse_2D/simulation2a/forward_model/jcm"
        phi=0
        job_ids = []

        # fix the keys.
        calc_keys = default_keys.copy()
        calc_keys.update(keys)

        # rotation scan.
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

        # wait for the calculations.
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

        # transform to intensity map.
        decayswitch=0
        Qxmat, Qzmat, Intensitymat = IntensityFEM2(thetaarray, config, results, decayswitch)
        matmeta=GetMatmeta(thetaarray, config, Qxmat, Qzmat, Intensitymat)
        
        index1=config.dimqxbranch+1
        index2=2*config.dimqxbranch+1
        #Qzdomain=matmeta[:,index1:index2,1]
        intmodel=matmeta[:,index1:index2,2]
        print("Forward model evaluated.")
        return intmodel




