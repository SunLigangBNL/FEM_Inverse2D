"""
This file contains all the previous work on CDSAXS by using JCMsuite.
"""
import sys
import os
import numpy as np
import time

jcm_root = "/home/sun2024/JCM_2024" # -> set JCMROOT installation directory.
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
jcmwave.info()

from forward_model.modules.kvectransform import kvectrans1a, kvectrans2, kvectrans2a
#from forward_model.modules.keys import keys as default_keys

class ForwardModel:
    def __init__(self, NThreads=1, Multiplicity=1):
        #self._data = load_data()
        # Start daemon at local machine
        jcmwave.daemon.shutdown()
        jcmwave.daemon.add_workstation(
            Hostname = "localhost",
            NThreads = NThreads,
            Multiplicity = Multiplicity)

    # Call the main FEM solver to solve the scattering problem.
    # input: thetaarray, folderinput, folderoutput, savenear, savefar.
    # output: solution in k space saved in txt format.
    #def EMScattering(self, thetaarray, dirinput, diroutput, savenear, savefar):
    def EMScattering(self, working_dir=None):
        #calc_keys = default_keys.copy()
        #calc_keys.update(keys)
        #totaltime1=time.time()

        job_ids=[]
        #os.makedirs(diroutput, exist_ok=True)
        thetaarray=[-30,0]
        phi=0
        keys={}

        for theta in thetaarray:
            print("*****************************************************************************")
            print("current thata : ",theta)
            #these_calc_keys = calc_keys.copy()
            #these_calc_keys["theta"] = theta
            #these_calc_keys["phi"] = phi
            keys['thetaphi']=np.array([theta, phi])

            if working_dir is not None:
                wd = os.path.join(working_dir,'phi{}_theta{}'.format(phi,theta))
            else:
                wd = None

            job_id = jcmwave.solve(
                "forward_model/jcm1/project.jcmpt",
                keys,
                temporary=(working_dir is None),
                working_dir=wd
            )
            job_ids.append(job_id)

         # wait for the calculations
        #results, logs = jcmwave.daemon.wait(job_ids=job_ids, verbose=False)
        results, logs = jcmwave.daemon.wait()
        return results, logs