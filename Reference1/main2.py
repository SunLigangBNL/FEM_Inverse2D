import sys
import os
jcm_root = "/home/sun2024/JCM_2024" # -> set JCMROOT installation directory
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
jcmwave.info()
import numpy as np
import time

# User interface area.
totaltime1=time.time()

thetaarray=np.linspace(-20,30,6)
print("thetaarray:",thetaarray)

keys={}
phi=0
folderoutput="Output_all"
os.makedirs(folderoutput, exist_ok=True)

for i in range(len(thetaarray)):
    theta=thetaarray[i]
    print("*****************************************************************************")
    print("current thata : ",theta)
    foldername = os.path.join(folderoutput, f"theta_{theta:.2f}")
    os.makedirs(foldername, exist_ok=True)
    
    keys['thetaphi']=np.array([theta, phi])
    temptime1=time.time()
    results=jcmwave.solve('project.jcmp', keys)
    temptime2=time.time()
    temptime = temptime2-temptime1
    print(f"Solution time in current iteration: {temptime:.2f} seconds")
    
    #nearresult=jcmwave.loadcartesianfields('project_results/nearfield.jcm')
    farresult = jcmwave.loadtable('project_results/fourier_modes_minus_y.jcm')
    
    #xarray=nearresult["X"]
    #yarray=nearresult["Y"]
    #zarray=nearresult["Z"]
    #fieldarray=nearresult["field"]
    
    #exspat=fieldarray[0][:,:,0]
    #eyspat=fieldarray[0][:,:,1]
    #ezspat=fieldarray[0][:,:,2]

    karray=farresult["K"]
    earray=farresult["ElectricFieldStrength"][0]
    
    # Here we only export the result data as what they are in JCMcontrol. 
    # Further processing is needed for correct coordinate transformation.

    kxarray=karray[:,0]
    kyarray=karray[:,1]
    kzarray=karray[:,2]

    exspec=earray[:,0]
    eyspec=earray[:,1]
    ezspec=earray[:,2]

    #np.savetxt(os.path.join(foldername, 'xarray.txt'), xarray, delimiter='\t', fmt='%.18e')
    #np.savetxt(os.path.join(foldername, 'yarray.txt'), yarray, delimiter='\t', fmt='%.18e')
    #np.savetxt(os.path.join(foldername, 'zarray.txt'), zarray, delimiter='\t', fmt='%.18e')

    #np.savetxt(os.path.join(foldername, 'exarray.txt'), exspat, delimiter='\t', fmt='%.18e')
    #np.savetxt(os.path.join(foldername, 'eyarray.txt'), eyspat, delimiter='\t', fmt='%.18e')
    #np.savetxt(os.path.join(foldername, 'ezarray.txt'), ezspat, delimiter='\t', fmt='%.18e')

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

