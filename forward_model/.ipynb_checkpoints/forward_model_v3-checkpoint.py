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
    def EMScattering(self, thetaarray, dirinput, diroutput, savenear, savefar, working_dir=None):

        totaltime1=time.time()
        job_ids=[]
        os.makedirs(diroutput, exist_ok=True)
        phi=0
        keys={}

        for theta in thetaarray:
            print("*********************************************")
            print("current thata : ",theta)
            keys['thetaphi']=np.array([theta, phi])

            if working_dir is not None:
                wd = os.path.join(working_dir,'phi{}_theta{}'.format(phi,theta))
            else:
                wd = None
                
            temptime1=time.time()
            job_id = jcmwave.solve(os.path.join(dirinput, "project.jcmp"),keys,temporary=(working_dir is None),working_dir=wd)    
            temptime2=time.time()
            temptime = temptime2-temptime1
            print(f"Solution time in current iteration: {temptime:.2f} seconds")            
            job_ids.append(job_id)

         # wait for the calculations
        results, logs = jcmwave.daemon.wait(job_ids=job_ids, verbose=False)

        # export the solution to txt files.
        for i in range(len(thetaarray)):
            theta=thetaarray[i]
            foldername = os.path.join(diroutput, f"theta_{theta:.2f}")
            os.makedirs(foldername, exist_ok=True)
            result=results[i]

            if savefar==True:
                farresult=result[3] # corresponding to the "fourier_modes_minus_y.jcm" in the project.jcmp file.
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

            if savenear == True:
                nearresult=result[1] # corresponding to the "nearfield.jcm" in the project.jcmp file.
            
                xarray=nearresult["X"]
                yarray=nearresult["Y"]
                zarray=nearresult["Z"]
                fieldarray=nearresult["field"]
                
                exspat=fieldarray[0][:,:,0]
                eyspat=fieldarray[0][:,:,1]
                ezspat=fieldarray[0][:,:,2]

                np.savetxt(os.path.join(foldername, 'xarray.txt'), xarray, delimiter='\t', fmt='%.18e')
                np.savetxt(os.path.join(foldername, 'yarray.txt'), yarray, delimiter='\t', fmt='%.18e')
                np.savetxt(os.path.join(foldername, 'zarray.txt'), zarray, delimiter='\t', fmt='%.18e')
            
                np.savetxt(os.path.join(foldername, 'exspat.txt'), exspat, delimiter='\t', fmt='%.18e')
                np.savetxt(os.path.join(foldername, 'eyspat.txt'), eyspat, delimiter='\t', fmt='%.18e')
                np.savetxt(os.path.join(foldername, 'ezspat.txt'), ezspat, delimiter='\t', fmt='%.18e')

            #log=np.array(logs[i])
            #np.savetxt(os.path.join(foldername, 'log.txt'), log, delimiter='\t', fmt='%.18e')

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

        totaltime2=time.time()
        totaltime = totaltime2 - totaltime1
        print(f"Total computation time: {totaltime:.2f} seconds")
        return results, logs






