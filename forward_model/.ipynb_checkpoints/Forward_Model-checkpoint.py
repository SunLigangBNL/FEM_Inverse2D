import os
import numpy as np
from contextlib import redirect_stdout
import time

import jcmwave
jcmwave.info()

from forward_model.modules.FEMProcessing import IntensityFEM, GetMatmeta
from forward_model.modules.BAProcessing import IntensityBA
from forward_model.modules.InterfaceTransform import karraytrans, FresnelFun
from forward_model.modules.keys import keys as default_keys
from forward_model.modules.utils import sortfunFEM

from pathlib import Path

class ForwardModel:
    def __init__(self, config, compenv, tempworkdir):
        self.config=config
        
        # Start daemon at local machine
        self.compenv=compenv
        jcmwave.daemon.shutdown()
        jcmwave.daemon.add_workstation(
            Hostname = "localhost",
            NThreads = compenv["NThreads"],
            Multiplicity = compenv["Multiplicity"])

        # # Start daemon on cluster.
        # self.compenv=compenv
        # jcmwave.daemon.shutdown()
        # queue_id = jcmwave.daemon.add_queue(**compenv)
        # print("Daemon Information after Add Queue:")
        # jcmwave.daemon.resource_info()
        
        # set up directories.
        BASE_DIR = Path(__file__).resolve().parent
        self.dir_jcm=BASE_DIR/"jcm"
        self.tempworkdir=tempworkdir
        os.environ['MODULES_DIR'] = str(BASE_DIR/"modules")
    
    # Compute reference: Part 1.
    # Input: keys, config, export directory.
    # Output: Output_all folder.
    def CompReference(self, keys, config, directory):

        totaltime1=time.time()

        # set up directories.
        folderoutput = os.path.join(directory, "Output_all")
        os.makedirs(folderoutput, exist_ok=True)
        
        # preparation.
        calc_keys=default_keys.copy()
        calc_keys.update(keys)
        print("current key under calculation:",calc_keys)
        
        job_ids=[]
        thetaarray=config.source.thetaarray
        dir_project_file=os.path.join(self.dir_jcm, "project.jcmpt")
        
        for theta in thetaarray:
            these_calc_keys=calc_keys.copy()
            these_calc_keys["theta"]=theta

            work_dir = os.path.join(self.tempworkdir, f"theta_{theta:.2f}")
            os.makedirs(work_dir, exist_ok=True)

            job_id=jcmwave.solve(dir_project_file, keys=these_calc_keys, temporary=False, working_dir=work_dir)
            job_ids.append(job_id)
        
        job_statuses = jcmwave.daemon.status(job_ids)
        #print("All status",job_statuses)

        results, logs = jcmwave.daemon.wait(job_ids = job_ids)

        totaltime2=time.time()
        totaltime = totaltime2-totaltime1
        print(f"Time of computing: {totaltime:.1f} seconds.")
        
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

        #print("*********************************** FEM Computing of Refernece Done ******************************************")


    # Compute reference: Part 2.
    # Input: Output_all folder.
    # Output: Qzmat, intrefmat, measurementmat, uncertaintymat.
    # Remark: the most interesting part of this function is how to generate the uncertainty of numerical reference.
    def CompReference2(self, config, directory):
        folderoutput=os.path.join(directory, "Output_all")
        Qxmat, Qzmat, Intensitymat, Psiradmat = IntensityFEM(config, folderoutput)
        matmeta = GetMatmeta(config, Qxmat, Qzmat, Intensitymat, Psiradmat)

        # export for BA.
        np.save(directory/"Qxmat.npy", Qxmat)
        np.save(directory/"Qzmat.npy", Qzmat)
        #np.save(directory/"Psiradmat.npy", Psiradmat)
        np.save(directory/"matmeta.npy", matmeta)
        
        qxindexarray=config.centerindex+config.qxpeakarray
        Qzmat2=matmeta[:,qxindexarray,1] # Processed Qz dataset.
        intrefmat=matmeta[:,qxindexarray,2] # True intensity without any noise. We will use it to generate the measurement data.

        # # Version 3: truncated Gaussian distribution.
        # mumat=np.zeros_like(intrefmat)
        # sigmafloor=(1e-6*np.max(intrefmat))**2   # manually introduce a floor for sigma: 1e-6 of max intensity.
        # afac=0.02 # more or less 1% of the reference data.
        # bfac=1e-16 # little positive background noise.
        # sigmamat=(afac*intrefmat)**2+bfac**2
        # sigmamat=np.maximum(sigmamat, sigmafloor)
        # muarray=mumat.ravel(order='C')
        # sigmaarray=sigmamat.ravel(order='C')
        # covmat=np.diag(sigmaarray)
        # noisearray=np.random.multivariate_normal(muarray, covmat)
        # noisemat=noisearray.reshape(mumat.shape,order='C')
        # measurementmat=intrefmat+noisemat # generate the measurement data containing noise.
        # measurementmat=np.maximum(measurementmat, 0.0) # to avoid negative values.
        # uncertaintymat=np.sqrt(sigmamat) # associated uncertainties.
        # return Qzmat, intrefmat, measurementmat, uncertaintymat

        # Version 4: Poisson distribution.
        gfac=1e7
        bgcounts=1e1
        lambdamat=gfac*intrefmat+bgcounts
        photoncounts=np.random.poisson(lambdamat)
        measurementmat=(photoncounts-bgcounts)/gfac
        uncertaintymat=np.sqrt(lambdamat)/gfac
        
        return Qzmat2, intrefmat, measurementmat, uncertaintymat        

    # Main forward model of FEM.
    # Input: single key, config, Qzgridmat.
    # Output: intensity matrix with standard dimension, after single point normalization.
    # Remark: DW decay was only applied to the observation peaks.
    def ModelEvaluateFEM(self, keys, config, Qzgridmat):

        thetaarray=config.source.thetaarray
        qxrefarray=np.linspace(-config.dimqxbranch*config.deltaqx,config.dimqxbranch*config.deltaqx,2*config.dimqxbranch+1)
        intmatFEM=np.zeros((len(thetaarray),len(qxrefarray)))

        # preparation.
        calc_keys=default_keys.copy()
        calc_keys.update(keys)
        #print("current key under calculation:",calc_keys)
        
        job_ids=[]
        dir_project_file=os.path.join(self.dir_jcm, "project.jcmpt")

        for theta in thetaarray:
            these_calc_keys=calc_keys.copy()
            these_calc_keys["theta"]=theta
            work_dir = os.path.join(self.tempworkdir, f"theta_{theta:.2f}")
            os.makedirs(work_dir, exist_ok=True)
            job_id=jcmwave.solve(dir_project_file, keys=these_calc_keys, temporary=False, working_dir=work_dir)
            job_ids.append(job_id)
        
        job_statuses = jcmwave.daemon.status(job_ids)
        print("All status",job_statuses)
        
        results, logs = jcmwave.daemon.wait(job_ids = job_ids)

        Qxmat, Qzmat, Intensitymat, Psiradmat = self.IntensityFEM2(config, results) #**************
        matmetaFEM = GetMatmeta(config, Qxmat, Qzmat, Intensitymat, Psiradmat)

        for i in range(len(qxrefarray)):
            matslice=matmetaFEM[:, i, :]
            Qzarray2, intarray2, thetaarray2, psiarray2=sortfunFEM(matslice)
            Qzgridmat[:,i]=Qzarray2
            intmatFEM[:,i]=intarray2

        # # Trun off the Angle-dependent factor.
        # print("ADF turned off.")
        # for i in range(len(qxrefarray)):
        #     if i==config.centerindex:
        #         continue
        #     Qxvalue=qxrefarray[i]
        #     Qzarray=Qzgridmat[:,i]
        #     Qlength=np.sqrt(Qxvalue**2+Qzarray**2)
        #     tempfac1=(Qlength/(2*config.source.k0))**2
        #     anglefac2=(1-2*tempfac1)**2/((Qxvalue/Qlength)**2-tempfac1)
        #     intmatFEM[:,i]=intmatFEM[:,i]/anglefac2

        # DW decay modification and scaling modification.
        if config.decayswitch==1:
            print("DW factors applied.")
            dwfacx2=keys['dwfacx']
            dwfacz2=keys['dwfacz']
            for i in range(len(qxrefarray)):
                if i==config.centerindex:
                    continue
                Qxvalue=qxrefarray[i]
                Qzarray=Qzgridmat[:,i]
                DWfacarray=np.exp(-((Qxvalue*1e-10*dwfacx2)**2+(Qzarray*1e-10*dwfacz2)**2))
                DWfacarray = np.clip(DWfacarray, 1e-300, None)  # avoid exact zeros
                if np.min(DWfacarray) < 1e-300:
                    print("min of DWfacarray (clipped):", np.min(DWfacarray))
                intmatFEM[:,i]=DWfacarray*intmatFEM[:,i]
            
        # normalization.
        if config.NLswitch==1:
            print("Intensity normalized.")
            Qzarraytemp=Qzgridmat[:,config.centerindex+1]
            Intarraytemp=intmatFEM[:,config.centerindex+1]
            mask2=~(Qzarraytemp==0)
            QzFEMarray=Qzarraytemp[mask2]
            IntFEMarray=Intarraytemp[mask2]
            normfacFEM=np.interp(0.0, QzFEMarray, IntFEMarray)
            intmatFEM=intmatFEM/normfacFEM

        print("************************** Forward Model FEM Evaluated *******************************")
        
        return intmatFEM

    # Main forward model of BA.
    # Input: single key, config, externally computed Qzmat (independent from FEM solver).
    # Output: intensity matrix with standard dimension, after single point normalization.
    # Remark: DW decay was only applied to the observation peaks.
    def ModelEvaluateBA(self, keys, config, Qzgridmat):
        
        #thetaarray=config.source.thetaarray
        qxrefarray=np.linspace(-config.dimqxbranch*config.deltaqx,config.dimqxbranch*config.deltaqx,2*config.dimqxbranch+1)
        intmatBA=IntensityBA(keys, config, Qzgridmat)     

        # DW decay modification and scaling modification.
        if config.decayswitch==1:
            print("DW factors applied.")
            dwfacx2=keys['dwfacx']
            dwfacz2=keys['dwfacz']
            for i in range(len(qxrefarray)):
                if i==config.centerindex:
                    continue
                Qxvalue=qxrefarray[i]
                Qzarray=Qzgridmat[:,i]
                DWfacarray=np.exp(-((Qxvalue*1e-10*dwfacx2)**2+(Qzarray*1e-10*dwfacz2)**2))
                DWfacarray = np.clip(DWfacarray, 1e-300, None)  # avoid exact zeros
                if np.min(DWfacarray) < 1e-300:
                    print("min of DWfacarray (clipped):", np.min(DWfacarray))
                intmatBA[:,i]=DWfacarray*intmatBA[:,i]
                
        
        # normalization.
        if config.NLswitch==1:
            print("Intensity normalized.")
            Qzarraytemp=Qzgridmat[:,config.centerindex-2] # based on the -2 peak.
            Intarraytemp=intmatBA[:,config.centerindex-2]
            mask2=~(Qzarraytemp==0)
            QzBAarray=Qzarraytemp[mask2]
            IntBAarray=Intarraytemp[mask2]
            normfacBA=np.interp(0.0, QzBAarray, IntBAarray)
            intmatBA=intmatBA/normfacBA
        
        print("*************************** Forward Model BA Evaluated **********************************")
        return intmatBA


    # Function from JCM results into matmeta directly, without exporting.
    # This function will only be called within ModelEvaluate function.
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
            #intarray=temparray1*temparray2/(config.source.E0**2*config.source.k0*np.abs(np.cos(thetarad)))
            intarray=temparray1*temparray2/(config.source.k0*np.abs(np.cos(thetarad)))
    
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
