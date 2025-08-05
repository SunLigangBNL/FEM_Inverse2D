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
from forward_model.modules.utils import selectsortfun2
from forward_model.modules.InterfaceTransform import karraytrans, FresnelFun
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
        print("current key under calculation:",calc_keys)
        
        phi=0

        # jcmwave.daemon.shutdown()
        # jcmwave.daemon.add_workstation(Hostname ='localhost', Multiplicity = 1, NThreads = 32)
        
        job_ids = []

        #currentdir = os.getcwd()
        tempdir = Path(__file__).resolve()
        currentdir = tempdir.parent
        project_pt = os.path.join(currentdir, "jcm", "project.jcmpt")
        
        for theta in thetaarray:
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

        # # checking area:
        # print("Current cd:", calc_keys["cd"])
        # print("Current h:", calc_keys["h"])

        totaltime2=time.time()
        totaltime = totaltime2 - totaltime1
        print("*****************************************************************************")
        print(f"Total computation time: {totaltime:.2f} seconds")


    # Compute reference: Part 2.
    # Input: Output_all folder.
    # Output: intmean, intuncertainty.
    def CompReference2(self, config, directory, qxpeakarray):
        folderoutput=directory
        
        Qxmat, Qzmat, Intensitymat, Psiradmat = IntensityFEM(config, folderoutput)
        matmeta = GetMatmeta(config, Qxmat, Qzmat, Intensitymat, Psiradmat)

        qxindexarray=config.centerindex+qxpeakarray
        Qzmat=matmeta[:,qxindexarray,1]
        intmeanmat=matmeta[:,qxindexarray,2]
        #intunctymat=intmeanmat*np.random.uniform(1,3)*0.0001+0.000001
        randmat=np.random.uniform(0.05, 0.2, size=intmeanmat.shape)
        intunctymat=intmeanmat*randmat+1e-12 # this is sigma!
        print("Numerical reference generated.")
        
        return Qzmat, intmeanmat, intunctymat
        


    # Main forward model.
    # Input: single key (fixed problem setup).
    # Output: Qzlist, intmeanlist, intunctylist
    def ModelEvaluate(self, keys, config, directory, qxpeakarray):
        os.makedirs(directory, exist_ok=True)

        folderworkdir=os.path.join(directory, "workdir_temp")
        os.makedirs(folderworkdir, exist_ok=True)

        # preparation.
        thetaarray=config.source.thetaarray
        dir_jcm=BASE_DIR/"jcm"

        calc_keys = default_keys.copy()
        calc_keys.update(keys)

        print("current key under calculation:",calc_keys)
        phi=0

        # jcmwave.daemon.shutdown()
        # jcmwave.daemon.add_workstation(Hostname ='localhost', Multiplicity = 1, NThreads = 32)
        
        job_ids = []

        tempdir = Path(__file__).resolve()
        currentdir = tempdir.parent
        project_pt = os.path.join(currentdir, "jcm", "project.jcmpt")
        
        for theta in thetaarray:
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
        
        print("Forward model evaluated.")

        Qxmat, Qzmat, Intensitymat, Psiradmat = self.IntensityFEM2(config, results) #**************
        matmeta = GetMatmeta(config, Qxmat, Qzmat, Intensitymat, Psiradmat)

        qxindexarray=config.centerindex+qxpeakarray
        intmeanmat=matmeta[:,qxindexarray,2]
            
        return intmeanmat

    # Function from JCM results into matmeta directly, without exporting.
    # This function will only be called with ModelEvaluate function.

    def IntensityFEM2(self, config, results):
        Qxmat=[]
        Qzmat=[]
        Intensitymat=[]
        Psiradmat=[]
        thetaarray=config.source.thetaarray
        for i in range(len(thetaarray)):
            theta=thetaarray[i]
            thetarad=np.radians(theta)
            # fill in karray and earray from results directly.
            farresult=results[i][3] #************************** careful!!
            karray=farresult["K"]
            earray=farresult["ElectricFieldStrength"][0]
    
            kxarray=np.real(karray[:,0])
            kyarray=np.real(karray[:,2])
            kzarray=np.real(karray[:,1]) #******
    
            exspec=earray[:,0]
            eyspec=earray[:,2]
            ezspec=earray[:,1] #******
    
            # The rest part is copied from the main IntensityFEM function.
            
            # Release diffracted waves from substrate into free space.
            # layerrefindexarray[-1] is always air.
            refindexlastlayer=np.real(config.material.layerrefindexarray[-2]) # Must take the real part to coinside with the JCM setting.
            refindexair=1.0
            kxarray3, kzarray3, diffradarray = karraytrans(kxarrayin=kxarray, kzarrayin=kzarray, refindexin=refindexlastlayer,
                                                           refindexout=refindexair, knormout=config.source.k0) 
    
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





