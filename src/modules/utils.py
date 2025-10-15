"""
This file is a library of the defined functions (small pieces).
"""

import numpy as np

# Coordinate transformation of the vertices of all scatterers.
# Input: [x,z] coordinates under JCM frame of all vertices in a list format.
# Output: a list of processed coorddinates of all vertices under sample frame.
# Remark: re-order the vertices as a preparation for Green's theorem.
def coordinate_transform(vertexmatlist1):
    vertexmatlist2=[]
    for vertexmat1 in vertexmatlist1:
        # step 1: transform to sample frame.
        vertexmat1a=np.column_stack((vertexmat1[:, 0], -vertexmat1[:, 1]))
        Nr=len(vertexmat1a) # number of vertices of the current scatterer.
        # step 2: reorder.
        tempz = vertexmat1a[:, 1].min()
        tempindex1 = np.where(vertexmat1a[:, 1] == tempz)[0]
        tempindex2 = tempindex1[np.argmin(vertexmat1a[tempindex1, 0])]
        vertexmat1b=np.array([vertexmat1a[(tempindex2-i)% Nr] for i in range(Nr)])
        vertexmatlist2.append(vertexmat1b)
    return vertexmatlist2
    
"""
Given mat_slice of shape (N, 5) for a fixed qx:
  mat_slice[:,1] = Qz
  mat_slice[:,2] = intensity
  mat_slice[:,3] = theta
  mat_slice[:,4] = psirad
Returns the filtered & sorted arrays (Qz, intensity, theta, psirad).
"""
def selectsortfun(matslice: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    Qzarray = matslice[:, 1]
    intarray = matslice[:, 2]
    thetaarray = matslice[:, 3]
    psiarray = matslice[:, 4]

    order = np.argsort(Qzarray)
    Qzarray2=Qzarray[order]
    intarray2=intarray[order]
    thetaarray2=thetaarray[order]
    psiarray2=psiarray[order]

    mask = ~((Qzarray2 == 0) & (intarray2 == 0))
    Qzarray3 = Qzarray2[mask]
    intarray3 = intarray2[mask]
    thetaarray3 = thetaarray2[mask]
    psiarray3 = psiarray2[mask]

    return Qzarray3, intarray3, thetaarray3, psiarray3

# A better version.
def selectsortfun2(matslice: np.ndarray):

    Qzarray = matslice[:, 1]
    intarray = matslice[:, 2]
    thetaarray = matslice[:, 3]
    psiarray = matslice[:, 4]

    order = np.argsort(Qzarray)
    Qzarray2=Qzarray[order]
    intarray2=intarray[order]
    thetaarray2=thetaarray[order]
    psiarray2=psiarray[order]

    # mask = ~((Qzarray2 == 0) & (intarray2 == 0))
    # Qzarray3 = Qzarray2[mask]
    # intarray3 = intarray2[mask]
    # thetaarray3 = thetaarray2[mask]
    # psiarray3 = psiarray2[mask]

    return Qzarray2, intarray2, thetaarray2, psiarray2


