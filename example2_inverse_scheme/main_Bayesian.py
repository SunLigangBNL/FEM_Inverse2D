# Goal: example to show how to use the inverse code.
# By: Ligang Sun.
# Date: 10/15/2025.

import os
import sys
import numpy as np
from pathlib import Path
import shutil
import time
import pickle
import uuid

# set up directory.
jcm_root = "/sdcc/u/lsun1/JCM_2025" #************* update.
jcm_optimizer_path = "/sdcc/u/lsun1/JCMoptimizer2" #*************
source_path="/sdcc/u/lsun1/Gitspace/FEM_Inverse2D" #*************

sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))

sys.path.insert(0, os.path.join(jcm_optimizer_path, "interface", "python"))
from jcmoptimizer import Server, Client, Study, Observation

sys.path.append(source_path)
from src.Forward_Model import ForwardModel
from src.modules import SourceConfig, GeometryConfig, Scatterer, MaterialConfig, CDSAXSConfig
from src.modules import IntensityFEM, GetMatmeta, selectsortfun, coordinate_transform

# Static library:
# rotation arrays:
thetaarray1=[-30,0]
thetaarray4=np.linspace(-60,60,121)
thetaarray5=np.linspace(-60,60,61)
thetaarray6=np.linspace(-60,60,31)
thetaarray7=np.linspace(-60,60,16)
thetaarray8=np.linspace(-60,60,7)

nanofac=1e-9

# material properties:
# wavelength: 1 nm.
epsr_Si=9.99435131564694e-1+3.278431066866658e-5j
epsr_W=9.97360870441372e-1+7.893104295072442e-4j

# define light source.
thetaarray=thetaarray7 #*************
lambda0=1e-9 # free space wavelength.*************
E0=1.0 # amplitude of incident electric field.
source=SourceConfig(thetaarray=thetaarray, lambda0=lambda0, E0=E0)

pitch=60*nanofac #*************
dimqxbranch=7 #*************
qxpeakarray=np.arange(0,7,1)
qxpeakarray=qxpeakarray[qxpeakarray != 0]
print("Qx peak order considered:",qxpeakarray)

# define layout.
# (1) background layered medium. Sample frame. Remark: layers are always seperated by different background materials.
layerzarray=np.array([0, 2e5])*nanofac # z-coordinates of layer interfaces.
scatterlayer=[0] # grating structure is on the top layer so index is 0.
scattercoordmat=[]
# (2) JCM frame: [x,z] corrdinates of the vertices of all scatterers. Copy from layout.jcm.
# Order: JCM layout top layer (high JCM z) first. From left to right.
vertexmat1 = np.array([[-15, 0], [15, 0],[12, 40],[-12, 40]])*nanofac
vertexmatlist=[vertexmat1] # JCM frame.
scattercoordmat.append(coordinate_transform(vertexmatlist)) # sample frame.
geometry=GeometryConfig(scattercoordmat=scattercoordmat, layerzarray=layerzarray, scatterlayer=scatterlayer)

# assign material property: JCM layout, from top to bottom, from left to right.
layerepsrarray=[1.0, epsr_Si, 1.0] # ****************
scatterepsrmat=[]
scatterepsrmat.append([epsr_W]) # ***************
material=MaterialConfig(layerepsrarray=layerepsrarray,scatterepsrmat=scatterepsrmat)

# generate the main config file.
config=CDSAXSConfig(pitch=pitch, dimqxbranch=dimqxbranch, source=source, geometry=geometry, material=material)

qxrefarray=np.linspace(-config.dimqxbranch*config.deltaqx,config.dimqxbranch*config.deltaqx,2*config.dimqxbranch+1)
qxindexarray=config.centerindex+qxpeakarray

# Part 1: Forward model initialization.
dir_source=Path().resolve()
dir_ref = dir_source/"Reference" #********
os.makedirs(dir_ref, exist_ok=True)
dir_model = dir_source/"Model_results" #**********
os.makedirs(dir_model, exist_ok=True)

dir_rltverr = os.path.join(dir_model, "converge_rltverr_report")
os.makedirs(dir_rltverr, exist_ok=True)

# computation environment.
compenv={
    'Hostname'      : 'localhost',
    'JCMROOT'       : jcm_root,
    'Type'          : 'Slurm',
    'Account'       : 'cfntest',
    'PartitionName' : 'cfn',
    'Time'          : 30, 
    'Multiplicity'  : 8, #*************
    'NNodes'        : 1,
    'NTasks'        : 2,
    'NTasksPerNode' : 2,
    'NThreads'      : 6 #*************
}

fm=ForwardModel(config, compenv)

# true parameters.
refkeys={'cd':27, 'h': 40, 'swa': 85.710846671181, 'r_top': 3.2, 'r_bot':5.5, 'pitch': 60, 
'epsilonr_scat_re':np.real(epsr_W), 'epsilonr_scat_im':np.imag(epsr_W)}
print("refkeys=",refkeys)

# # Part 2: generate numerical reference (with Poisson noise term).
# ref=fm.CompReference(keys=refkeys,config=config,directory=dir_ref)
# Qzref, inttrueref, intmearef, intunctyref=fm.CompReference2(config=config,directory=dir_ref, qxpeakarray=qxpeakarray)

# # Part 2a: export the random measurement data to external file.
# np.save(dir_ref/"Qzref.npy", Qzref)
# np.save(dir_ref/"inttrueref.npy", inttrueref)
# np.save(dir_ref/"intmearef.npy", intmearef)
# np.save(dir_ref/"intunctyref.npy",intunctyref)

# Part 2b: import the previously generated measurement data and uncertainty.
intmearef=np.load(dir_ref/"intmearef.npy")
intunctyref=np.load(dir_ref/"intunctyref.npy")

# Part 3: prepare for the main Bayesian engine.
server = Server()
client = Client(server.host)

# define search space.
varyingrange=0.2
facarray1=np.array([0.81,0.88,0.94,0.97,0.85]) # for 0.2 varying.
print("lower bound of parameter space : ", facarray1)
facarray2=facarray1+varyingrange
print("upper bound of parameter space : ", facarray2)

design_space = [
{'name': 'cd', 'domain': [refkeys['cd']*facarray1[0], refkeys['cd']*facarray2[0]]},
{'name': 'h', 'domain': [refkeys['h']*facarray1[1], refkeys['h']*facarray2[1]]},
{'name': 'swa', 'domain': [refkeys['swa']*facarray1[2], refkeys['swa']*facarray2[2]]},
{'name': 'r_top', 'domain': [refkeys['r_top']*facarray1[3], refkeys['r_top']*facarray2[3]]},
{'name': 'r_bot', 'domain': [refkeys['r_bot']*facarray1[4], refkeys['r_bot']*facarray2[4]]}
]

param_names = {param['name'] : '${' + '}_{\\rm '.join(param['name'].split('_')) + '}$' for param in design_space}

# define fixed parameter.
environment = [
{'name': 'pitch', 'type': 'fixed', 'domain': refkeys['pitch']},
{'name': 'epsilonr_scat_re', 'type': 'fixed', 'domain': refkeys['epsilonr_scat_re']},
{'name': 'epsilonr_scat_im', 'type': 'fixed', 'domain': refkeys['epsilonr_scat_im']}
]

constraints = [
{'name': 'swa', 'expression': 'swa <= 90'},
{'name': 'r_top', 'expression': 'r_top + h/2/tan(PI*swa/180) <= cd/2'},
{'name': 'r_bot', 'expression': 'r_bot + cd/2 + h/2/tan(PI*swa/180) <= pitch/2'}
]

# every time we should use a new study_id.
study = client.create_study(
    design_space=design_space,
    environment=environment,
    constraints=constraints,
    driver="BayesianLeastSquares",
    name="GeoRecon",
    study_id=str(uuid.uuid4().hex[:4]),
    save_dir=os.getcwd()
)

# defined iteration parameters.
itnrmax1=60 # itnr limit to find opt. solution.
itnrmax2=600 # itnr limit to warm up the MCMC engine.
accgoal=1.0e-3
autostop=True

# define target_vector and uncertainty vector.
target_vector=[*[value for i in range(len(qxindexarray)) for value in intmearef[:, i]]]
uncertainty_vector=[*[value for i in range(len(qxindexarray)) for value in intunctyref[:, i]]]
refnorm=np.linalg.norm(np.array(target_vector))

study.configure(
    max_iter=itnrmax1,
    target_vector=target_vector,
    uncertainty_vector=uncertainty_vector
)

def evaluate(study, **kwargs):
    if not hasattr(evaluate, "counter"):
        evaluate.counter = 0
    evaluate.counter += 1

    observation=study.new_observation()
    intmodmat=fm.ModelEvaluate(kwargs, config=config, directory=dir_model, qxpeakarray=qxpeakarray)

    observation_vector=[*[value for i in range(len(qxpeakarray)) for value in intmodmat[:, i]]]
    observation.add(observation_vector)

    print("Relative error report:")
    rltverrarray=[]
    for i in range(len(qxpeakarray)):
        array1=intmearef[:,i]
        array2=intmodmat[:,i]
        rltverr=np.linalg.norm(array1-array2)/np.linalg.norm(array1)
        rltverrarray.append(rltverr)

    rltverrall2=np.mean(rltverrarray)
    print(f"Iteration {evaluate.counter}, average relative error = {rltverrall2:.4f}")

    fname=os.path.join(dir_rltverr, f"rltverr_run_{evaluate.counter}.txt")
    np.savetxt(fname, rltverrarray, fmt="%.6f")
    
    # Stop the engine when the threshold is reached.
    if autostop==True:
        if rltverrall2 < accgoal :
            print("****************************** Accuracy Goal Reached. Bayesian Engine Stopped. ********************************")
            study.configure(max_iter=evaluate.counter)
    return observation

study.set_evaluator(evaluate)


# Part 4: run the main Bayesian engine.
totaltime1=time.time()
study.run()
totaltime2=time.time()
totaltime = totaltime2-totaltime1
print("***************** Convergence Part Done !! ********************")
print(f"Solution time of convergence: {totaltime:.2f} seconds")

# analysis of the opt. solution.
best_sample = study.driver.best_sample
min_chisq = study.driver.min_objective
uncertainties = study.driver.uncertainties
print("Optimized key:",best_sample)
for param in design_space:
    name = param['name']
    print(f"  {name} = {best_sample[name]:.5f} +/- {uncertainties[name]:.5f}")
print(f"Reconstructed parameters with chi-squared value {min_chisq:.4e}:")

optkeys=dict(refkeys)
for k, v in best_sample.items():
    if k in optkeys:
        optkeys[k] = v
print("optkeys=",optkeys)

optuncty=dict(refkeys)
for k, v in uncertainties.items():
    if k in optkeys:
        optuncty[k] = v
print("optuncty=",optuncty)

np.save(dir_model/"optkeys.npy",optkeys)
np.save(dir_model/"optuncty.npy",optuncty)


# Part 5: MCMC analysis.

print("********************************** Warm up the MCMC Engine *********************************************")

totaltime1=time.time()

# From here we attack the MCMC procedure.
# Before running a Markov-chain Monte-Carlo (MCMC) sampling we converge the surrogate
# models by sampling around the minimum. To make the study more explorative, the
# scaling parameter is increased and the effective degrees of freedom is set to one.
autostop=False
study.configure(
    scaling=10.0,
    #effective_DOF=1.0,
    min_uncertainty=max(min_chisq*1e-8, 1e-8),
    max_iter=itnrmax2,
    min_val=0.0
)
study.run()

totaltime2=time.time()
totaltime = totaltime2-totaltime1

print(f"Burn-in phase of MCMC: {totaltime:.2f} seconds")

temp_dir=dir_model/"workdir_temp"
shutil.rmtree(temp_dir)


print("********************************** Start MCMC *********************************************")

# Run the MCMC sampling with 32 walkers
num_walkers, max_iter = 64, 50000 #************************
mcmc_result=study.driver.run_mcmc(
    rel_error=0.01,
    num_walkers=num_walkers,
    max_iter=max_iter,
    max_sigma_dist=3.0
)

outfile = dir_model/"mcmc_result.pkl"
with open(outfile, "wb") as f:
    pickle.dump(mcmc_result, f)

print(f"MCMC result saved to {outfile}")

print("********************************** Done !!! *********************************************")
