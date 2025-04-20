"""
This file contains all the defined parameters to set up the CDSAXS problem.
"""

import numpy as np

class ConfigParameter:
    def __init__(self, thetaarray, pitch, line, height, thickness, rr, epsilonr1, epsilonr2, lambda0, dimqxbranch):
        
        self.thetaarray=thetaarray
        
        self.pitch=pitch
        self.line=line
        self.height=height
        self.thickness=thickness
        self.rr=rr

        self.epsilonr1=epsilonr1 # grating line.
        self.epsilonr2=epsilonr2 # substrate.
        
        self.refindexi=1.0 # free space.
        self.refindext=np.sqrt(epsilonr2)
        
        self.lambda0=lambda0
        self.E0=1.0 # amplitude of the incident electric field.
        self.k0=2*np.pi/lambda0
        self.lambda1=lambda0/self.refindext
        self.k1=2*np.pi/self.lambda1
        self.dimqxbranch=dimqxbranch
