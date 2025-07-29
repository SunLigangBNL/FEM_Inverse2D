import sys
import os
import numpy as np
from contextlib import redirect_stdout
import time

#jcm_root = "/sdcc/u/lsun1/JCM_2025"
jcm_root = "/home/sun2024/JCM_2024"
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
#jcmwave.info()

from forward_model.modules.FEMProcessing import IntensityFEM, GetMatmeta
from forward_model.modules.utils import selectsortfun
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
        
    # Compute reference: Part 1.
    # Input: keys, config, export directory.
    # Output: Output_all folder.
    def CompReference(self, keys, config, directory):
        
        totaltime1=time.time()
        os.makedirs(directory, exist_ok=True)
        
        folderoutput = os.path.join(directory, "Output_all")
        os.makedirs(folderoutput, exist_ok=True)

        folderworkdir=os.path.join(directory, "workdir_temp")
        os.makedirs(folderworkdir, exist_ok=True)

        # preparation.
        thetaarray=config.source.thetaarray
        dir_jcm=BASE_DIR/"jcm"

        calc_keys = default_keys.copy()
        calc_keys.update(keys)
        phi=0

        jcmwave.daemon.shutdown()
        jcmwave.daemon.add_workstation(Hostname ='localhost', Multiplicity = 1, NThreads = 32)
        
        job_ids = []

        #currentdir = os.getcwd()
        tempdir = Path(__file__).resolve()
        currentdir = tempdir.parent
        project_pt = os.path.join(currentdir, "jcm", "project.jcmpt")
        
        for theta in thetaarray:
            if theta%10==0:
                print("current thata : ",theta)
            these_calc_keys = calc_keys.copy()
            these_calc_keys["theta"] = theta
            these_calc_keys["phi"] = phi
            work_dir = os.path.join(folderworkdir, f"theta_{theta:.2f}")
            os.makedirs(work_dir, exist_ok=True)
            job_id=jcmwave.solve(project_pt, keys=these_calc_keys, temporary=False, working_dir=work_dir)
            job_ids.append(job_id)
        
        job_statuses = jcmwave.daemon.status(job_ids)
        print("All status",job_statuses)

        results, logs = jcmwave.daemon.wait(job_ids = job_ids)

        print("Main computation done.")
        
        log_filename = os.path.join(folderoutput, "logs_all.txt")
        with open(log_filename, "w") as log_file, redirect_stdout(log_file):
            for i, log in enumerate(logs):
                print(f"Log {i}:")
                print(log["Log"]["Out"])
                print("*" * 100)

        print(f"All logs saved in {log_filename}.")

        # export solution results.
        for i in range(len(thetaarray)):
            theta=thetaarray[i]
            
            foldername = os.path.join(folderoutput, f"theta_{theta:.2f}")
            os.makedirs(foldername, exist_ok=True)
            
            farresult=results[i][3] #**************************
            
            karray=farresult["K"]
            earray=farresult["ElectricFieldStrength"][0]
        
            kxarray=karray[:,0]
            kyarray=karray[:,1]
            kzarray=karray[:,2]
        
            exspec=earray[:,0]
            eyspec=earray[:,1]
            ezspec=earray[:,2]
        
            np.savetxt(os.path.join(foldername, 'kxarray.txt'), kxarray, delimiter='\t', fmt='%.18e')
            np.savetxt(os.path.join(foldername, 'kyarray.txt'), kyarray, delimiter='\t', fmt='%.18e')
            np.savetxt(os.path.join(foldername, 'kzarray.txt'), kzarray, delimiter='\t', fmt='%.18e')
        
            np.savetxt(os.path.join(foldername, 'exspec.txt'), exspec, delimiter='\t', fmt='%.18e')
            np.savetxt(os.path.join(foldername, 'eyspec.txt'), eyspec, delimiter='\t', fmt='%.18e')
            np.savetxt(os.path.join(foldername, 'ezspec.txt'), ezspec, delimiter='\t', fmt='%.18e')

        print("Solution files exported.")

        totaltime2=time.time()
        totaltime = totaltime2 - totaltime1
        print("*****************************************************************************")
        print(f"Total computation time: {totaltime:.2f} seconds")


    # Compute reference: Part 2.
    # Input: Output_all folder.
    # Output: intmean, intuncertainty.
    def CompReference2(self, config, directory, qxpeakarray):
        folderoutput=directory
        qxindexarray=config.centerindex+qxpeakarray
        
        Qxmat, Qzmat, Intensitymat, Psiradmat = IntensityFEM(config, folderoutput)
        matmeta = GetMatmeta(config, Qxmat, Qzmat, Intensitymat, Psiradmat)

        Qzlist=[]
        intmeanlist=[]
        intunctylist=[]
        
        for index in qxindexarray:
            matslice = matmeta[:, index, :]
            Qzarray3, IntFEMarray, thetaarray3, psiradarray3 = selectsortfun(matslice)
            intunctyarray=IntFEMarray*np.random.uniform(0.1,0.4)+0.00001
            Qzlist.append(Qzarray3)
            intmeanlist.append(IntFEMarray)
            intunctylist.append(intunctyarray)

        return Qzlist, intmeanlist, intunctylist


    # # Main forward model.
    # # Input: single key (fixed problem setup).
    # # Process: rotation scan of multiple simulations.
    # # Output: corresponding intensity map.
    # def ModelEvaluate(self, keys, config, working_dir=None):

    #     # preparation.
    #     thetaarray=config.source.thetaarray
    #     directory=BASE_DIR/"jcm"
        
    #     phi=0
    #     job_ids = []

    #     # fix the keys.
    #     calc_keys = default_keys.copy()
    #     calc_keys.update(keys)

    #     # main rotation scan.
    #     for theta in thetaarray:
    #         if theta%10==0:
    #             print("current thata : ",theta)
    #         these_calc_keys = calc_keys.copy()
    #         these_calc_keys["theta"] = theta
    #         these_calc_keys["phi"] = phi

    #         if working_dir is not None:
    #             wd = os.path.join(working_dir,'phi{}_theta{}'.format(phi,theta))
    #         else:
    #             wd = None
                
    #         job_id = jcmwave.solve(os.path.join(directory, "project.jcmpt"),
    #                                keys=these_calc_keys,
    #                                temporary=(working_dir is None),
    #                                working_dir=wd)
    #         job_ids.append(job_id)

    #     # wait for the calculations.
    #     job_statuses = jcmwave.daemon.status(job_ids)
    #     print("All status",job_statuses)
    #     results, logs = jcmwave.daemon.wait(job_ids=job_ids, verbose=False)

    #     # export log files.
    #     log_filename = os.path.join(directory, "logs_all.txt")
    #     with open(log_filename, "w") as log_file, redirect_stdout(log_file):
    #         # anything printed in here goes to log_file only
    #         for i, log in enumerate(logs):
    #             print(f"Log {i}:")
    #             print(log["Log"]["Out"])
    #             print("*" * 100)

    #     print(f"All logs saved in {log_filename}.")

        

    #     # transform to intensity map.
    #     decayswitch=0
    #     Qxmat, Qzmat, Intensitymat = IntensityFEM2(config, results, decayswitch)
    #     matmeta=GetMatmeta(config, Qxmat, Qzmat, Intensitymat)
        
    #     index1=config.dimqxbranch+1
    #     index2=2*config.dimqxbranch+1
    #     #Qzdomain=matmeta[:,index1:index2,1]
    #     intmodel=matmeta[:,index1:index2,2]
    #     print("Forward model evaluated.")
    #     return intmodel




