import sys
import os
jcm_root = "/home/sun2024/JCM_2024" # update.
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
jcmwave.info()
import numpy as np
import time
from contextlib import redirect_stdout

# User interface area.
totaltime1=time.time()

thetaarray=np.linspace(-60,50,111)
print("thetaarray:",thetaarray)

keys={}
phi=0
folderoutput="Output_all"
os.makedirs(folderoutput, exist_ok=True)

jcmwave.daemon.shutdown()
jcmwave.daemon.add_workstation(Hostname ='localhost', Multiplicity = 1, NThreads = 32)

job_ids = []

for theta in thetaarray:
    keys['thetaphi']=np.array([theta, phi])
    job_id = jcmwave.solve('project.jcmp', keys=keys, temporary = True)
    job_ids.append(job_id)

job_statuses = jcmwave.daemon.status(job_ids)
print("All status",job_statuses)

results, logs = jcmwave.daemon.wait(job_ids = job_ids)

log_filename = os.path.join(folderoutput, "logs_all.txt")

with open(log_filename, "w") as log_file, redirect_stdout(log_file):
    for i, log in enumerate(logs):
        print(f"Log {i}:")
        print(log["Log"]["Out"])
        print("*" * 100)

print(f"All logs saved in {log_filename}.")

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

totaltime2=time.time()
totaltime = totaltime2 - totaltime1
print("*****************************************************************************")
print(f"Total computation time: {totaltime:.2f} seconds")


