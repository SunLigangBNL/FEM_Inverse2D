import sys
import os
import numpy as np
from contextlib import redirect_stdout

jcm_root = "/home/sun2024/JCM_2024"
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
#jcmwave.info()

from forward_model.modules.FEMProcessing import IntensityFEM, GetMatmeta
from forward_model.modules.keys import keys as default_keys


from pathlib import Path
BASE_DIR = Path(__file__).resolve().parent

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
    def CompRef(self, config, directory):
        decayswitch=0
        Qxmat, Qzmat, Intensitymat = IntensityFEM(config, directory, decayswitch)
        matmeta=GetMatmeta(config, Qxmat, Qzmat, Intensitymat)

        # select from the first peak.
        index1=config.dimqxbranch+1
        index2=2*config.dimqxbranch+1
        Qzdomain=matmeta[:,index1:index2,1]
        intmean=matmeta[:,index1:index2,2]
        intuncertainty=intmean*np.random.uniform(0.1,0.4)+0.00001
        print("Numerical reference generated.")
        return Qzdomain, intmean, intuncertainty

    # Main forward model.
    # Input: single key (fixed problem setup).
    # Process: rotation scan of multiple simulations.
    # Output: corresponding intensity map.
    def ModelEvaluate(self, keys, config, working_dir=None):

        # preparation.
        thetaarray=config.thetaarray
        directory=BASE_DIR/"jcm"
        
        phi=0
        job_ids = []

        # fix the keys.
        calc_keys = default_keys.copy()
        calc_keys.update(keys)

        # main rotation scan.
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
                
            job_id = jcmwave.solve(os.path.join(directory, "project.jcmpt"),
                                   keys=these_calc_keys,
                                   temporary=(working_dir is None),
                                   working_dir=wd)
            job_ids.append(job_id)

        # wait for the calculations.
        job_statuses = jcmwave.daemon.status(job_ids)
        print("All status",job_statuses)
        results, logs = jcmwave.daemon.wait(job_ids=job_ids, verbose=False)

        # # export the log file.
        # log_filename = os.path.join(directory, "logs_all.txt")
        # with open(log_filename, "w") as log_file:
        #     sys.stdout = log_file
        #     for i, log in enumerate(logs):
        #         print(f"Log {i}:")
        #         print(log["Log"]["Out"])
        #         print("*" * 100)
        #     sys.stdout = sys.__stdout__
        # print(f"All logs saved in {log_filename}.")

        # export the log file.
        # instead of manually doing sys.stdout = log_file, use:
        log_filename = os.path.join(directory, "logs_all.txt")
        with open(log_filename, "w") as log_file, redirect_stdout(log_file):
            # anything printed in here goes to log_file only
            for i, log in enumerate(logs):
                print(f"Log {i}:")
                print(log["Log"]["Out"])
                print("*" * 100)

        print(f"All logs saved in {log_filename}.")

        

        # transform to intensity map.
        decayswitch=0
        Qxmat, Qzmat, Intensitymat = IntensityFEM2(config, results, decayswitch)
        matmeta=GetMatmeta(config, Qxmat, Qzmat, Intensitymat)
        
        index1=config.dimqxbranch+1
        index2=2*config.dimqxbranch+1
        #Qzdomain=matmeta[:,index1:index2,1]
        intmodel=matmeta[:,index1:index2,2]
        print("Forward model evaluated.")
        return intmodel




