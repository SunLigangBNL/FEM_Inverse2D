"""
This file contains all the previous work on CDSAXS by using JCMsuite.
"""
import sys
import os
import numpy as np
import itertools
import time

jcm_root = "/home/sun2024/JCM_2024" # -> set JCMROOT installation directory.
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
jcmwave.info()

from forward_model.modules.kvectransform import kvectrans1a, kvectrans2, kvectrans2a
from forward_model.modules.keys import keys as default_keys
from forward_model.modules.load_data import load_data

class ConfigParameter:
    def __init__(self, pitch, line, height, relpermittivity, refindexi, refindext, lambda0, k0, lambda1, k1):
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


    def __repr__(self):
        params = ", ".join(f"{key}={value}" for key, value in self.__dict__.items())
        return f"ConfigParameter({params})"


class ForwardModel:
    def __init__(self, config, NThreads=1, Multiplicity=1):
        self._data = load_data()
        self.config=config
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
        print("Returning from EMScattering...")
        return results, logs

    # Compute far field intensity in the Q space.
    # Input: thetaarray.
    # Output: Qxmat, Qzmat, Intensitymat, matmeta.
    def SolTransform(self, thetaarray, config, diroutput, dimqxtemp):
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
            #if i%20==0:
            print("current thata : ",theta)
                
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

            # actually we do not need this filtering step.
            #indexlist2= [i for i, Qxvalue in enumerate(Qxarray3) if -ansbound<=Qxvalue<=ansbound]
            #Qxarray4= [Qxarray3[i] for i in indexlist2]
            #Qzarray4= [Qzarray3[i] for i in indexlist2]
            #intarray4=[intarray3[i] for i in indexlist2]
            
            Qxmat.append(Qxarray3)
            Qzmat.append(Qzarray3)
            Intensitymat.append(intarray3)
        
        Qxmat=np.array(Qxmat, dtype=object)
        Qzmat=np.array(Qzmat, dtype=object)
        Intensitymat=np.array(Intensitymat, dtype=object)
        #Intensitymat2=[np.log10(np.array(sublist, dtype=np.float64)) for sublist in Intensitymat]
        print("Intensity in Q-space computed.")
        
        periodicityx=config.pitch*10**10 # unit: angstrom.
        deltaqx=2*np.pi/periodicityx # unit: per angstrom.
        dimqx=2*dimqxtemp+1 # number of all interested Qx.
        print("number of interested Qx values in total: ", dimqx)
        qxrefarray=np.linspace(-dimqxtemp*deltaqx,dimqxtemp*deltaqx,dimqx)
        
        matmeta=np.zeros((len(thetaarray),dimqx,3)) # per angle, per Qx location, store (Qx, Qz, Intensity) values.
        for i in range(len(thetaarray)):
            if i%20==0:
                print("i=",i)
            qxarray=np.array(Qxmat[i])
            for j,qxvalue in enumerate(qxrefarray):
                index1 = np.argmin(np.abs(qxarray-qxvalue))
                if np.abs(qxarray[index1]-qxvalue)<0.5*deltaqx:
                    #print(qxarray[index1])
                    #print(qxvalue)
                    matmeta[i,j,0]=Qxmat[i][index1]
                    matmeta[i,j,1]=Qzmat[i][index1]
                    matmeta[i,j,2]=Intensitymat[i][index1]
        
        print("matmeta completed.")
        
        return Qxmat, Qzmat, Intensitymat, matmeta

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
            #if i%20==0:
            print("current thata : ",theta)
                
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
            if i%20==0:
                print("i=",i)
            qxarray=np.array(Qxmat[i])
            for j,qxvalue in enumerate(qxrefarray):
                index1 = np.argmin(np.abs(qxarray-qxvalue))
                if np.abs(qxarray[index1]-qxvalue)<0.5*deltaqx:
                    matmeta[i,j,0]=Qxmat[i][index1]
                    matmeta[i,j,1]=Qzmat[i][index1]
                    matmeta[i,j,2]=Intensitymat[i][index1]
        
        print("matmeta completed.")
        
        Qzmat=matmeta[:,:,1]
        int_mean=matmeta[:,:,2]
        int_uncertainty=int_mean*np.random.uniform(-0.02,0.02)
        
        return Qzmat, int_mean, int_uncertainty

    #def ModelEvaluate(self, thetaarray, config, diroutput, dimqxtemp):
    #    return result
    
    def ModelEvaluate(self, keys, dirinput, working_dir=None):
        """Run the numerical model and return the simulation results in a dictionary.
         
        :param keys: Argument keys for the simulation.
        :param working_dir: Directory, where simulation results are stored.
                     (default: temporary directory is used)

        """

        calc_keys = default_keys.copy()
        calc_keys.update(keys)
        
        import_results = load_data()

        thetas = self._data[0]["thetas"]
        thetas=thetas[2:5]
        phis = set(r["phi"] for r in self._data)

        t_p_combinations = list(itertools.product(thetas, phis))
        job_ids = []

        for i, (theta, phi) in enumerate(t_p_combinations):
            these_calc_keys = calc_keys.copy()
            these_calc_keys["theta"] = theta
            these_calc_keys["phi"] = phi
            print("theta=",theta)
            print("phi=",phi)
            if working_dir is not None:
                wd = os.path.join(working_dir,'phi{}_theta{}'.format(phi,theta))
            else:
                wd = None
            job_id = jcmwave.solve(os.path.join(dirinput, "project.jcmpt"),keys=these_calc_keys,temporary=(working_dir is None),working_dir=wd)
            job_ids.append(job_id)

        # wait for the calculations
        job_statuses = jcmwave.daemon.status(job_ids)
        print("All status",job_statuses)
        results, logs = jcmwave.daemon.wait()

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
























    

    def model_eval(self, keys, dirinput, working_dir=None):
        """Run the numerical model and return the simulation results in a dictionary.
         
        :param keys: Argument keys for the simulation.
        :param working_dir: Directory, where simulation results are stored.
                     (default: temporary directory is used)

        """

        calc_keys = default_keys.copy()
        calc_keys.update(keys)
        
        import_results = load_data()

        thetas = self._data[0]["thetas"]
        thetas=thetas[2:5]
        phis = set(r["phi"] for r in self._data)
        #print("*********")
        t_p_combinations = list(itertools.product(thetas, phis))

        job_ids = []

        #kk=0
        for i, (theta, phi) in enumerate(t_p_combinations):
            these_calc_keys = calc_keys.copy()
            these_calc_keys["theta"] = theta
            these_calc_keys["phi"] = phi
            print("theta=",theta)
            print("phi=",phi)
            if working_dir is not None:
                wd = os.path.join(working_dir,'phi{}_theta{}'.format(phi,theta))
            else:
                wd = None
            #job_id = jcmwave.solve("forward_model/jcm3/project.jcmpt",keys=these_calc_keys,temporary=(working_dir is None),working_dir=wd)
            job_id = jcmwave.solve(os.path.join(dirinput, "project.jcmpt"),keys=these_calc_keys,temporary=(working_dir is None),working_dir=wd)
            job_ids.append(job_id)

        # wait for the calculations
        #results, logs = jcmwave.daemon.wait(job_ids=job_ids, verbose=False)
        job_statuses = jcmwave.daemon.status(job_ids)
        print("All status",job_statuses)
        results, logs = jcmwave.daemon.wait()

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
        
        return results





















