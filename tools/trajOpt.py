import numpy as np
from orbits.orbit import Orbit
from tools.utils import uLambert

#FIXME: Seems to return low values 
def getC3(body, r2, tof):
    """
    Determines depart energy needed between a body and a destination
    INPUT :
        body    : The primary body we're departing from
        r2      : The destination (x,y,z) list
        tof     : The time of flight in seconds

    OUTPUT :
        c3      : The depart energy required
    """
    try:
        r1 = body.Orbit.cartesianCoordinates()
        v1 = body.Orbit.getSpeed()
        _,vE,_,_ = uLambert(r1, r2, tof, body.primBody)
        vE = np.linalg.norm(vE)
        print(vE)
        print(v1)
    except:
        return np.nan
    return (vE-v1)**2

#
# def getDepartEnergy(body1, body2, tof):
#     if body1.primBody != body2.primBody:
#         raise ValueError("The two bodies need the same primary body !")
#
#     body2.Orbit.updTime(tof)
#
#     v1 = body1.Orbit.getSpeed()
#     try:
#         r1 = np.array(body1.Orbit.cartesianCoordinates())
#         r2 = np.array(body2.Orbit.cartesianCoordinates())
#         _,vTransfer,_,_ = uLambert(r1, r2, tof, body1.primBody)
#         vTransfer = np.linalg.norm(vTransfer)
#     except:
#         body2.Orbit.updTime(-tof)
#         return np.nan
#
#     body2.Orbit.updTime(-tof)
#     return (vTransfer-v1)**2

def getPorkChop(body1, body2, departInterval, arrivalInterval, ITER):
    depArray, stepD = np.linspace(departInterval[0], departInterval[1], num = ITER, retstep=True)
    arrArray, stepA = np.linspace(arrivalInterval[1], arrivalInterval[0],num = ITER, retstep=True)

    retGrid = list()
    tofGrid = list()
    i=0

    for arr in np.nditer(arrArray):
        i+=1

        body1.Orbit.updTime(stepD)
        body2.Orbit.updTime(stepD)

        C3tof = list()
        tofList = list()

        for dep in np.nditer(depArray):
            tof = arr-dep
            if tof < 100*3600*24:
                C3 = np.nan
                C3tof.append(C3)
                tofList.append(tof/(3600*24))
                continue

            C3 = getDepartEnergy(body1, body2, tof)

            # if C3 is not np.nan:
            #     C3 = np.round(C3)

            C3tof.append(C3)
            tofList.append(tof/(3600*24))

        retGrid.append(C3tof)
        tofGrid.append(tofList)

        if i%5==0:
            print(str((i/ITER)*100)+"%")

    return np.array(retGrid, dtype=float), np.array(tofGrid)
