import numpy as np
from scipy.optimize import fmin

def sphereToCartesian2D(angle, distance):
    """
    INPUT :
        angle     - between (0, 2pi)
        distance  - between (0, inf)
    OUTPUT :
        tuple (x,y), 2d coordinates
    """
    x,y = 0,0
    x = distance*np.cos(angle)
    y = distance*np.sin(angle)

    return (x,y)

def sphereToCartesian3D(theta, phi, distance):
    """
    INPUT :
        theta       - between (0,2pi)   - azimuthal angle   (x-y plane from x  )
        phi         - between (0,pi)    - polar angle       (from z+ to direct.)
        distance    - between (0, inf)  - distance to origin
    OUTPUT :
        tuple (x,y,z), 3d coordinates
    """
    x = distance*np.cos(theta)*np.sin(phi)
    y = distance*np.sin(theta)*np.sin(phi)
    z = distance*np.cos(phi)
    return (x,y,z)

def meanToTrue(Me, ecc):
    """
    INPUT :
        Me          - between(0, 2pi)   - Mean anomaly
        ecc         - between(0, 1)     - eccentricity of the orbit
    OUTPUT :
        trueAn      - between(0 , 2pi)  - True anomaly


    Solves Kepler's equation for the true anomaly using Alg. 3.1

    """
    Me = Me%(2*np.pi)
    if Me==0:
        return 0

    if Me<np.pi:
        linEcc = Me+ecc/2
    else:
        linEcc = Me-ecc/2

    func = linEcc-ecc*np.sin(linEcc)-Me
    deriv = 1-ecc*np.cos(linEcc)
    ratio = func/deriv

    while np.abs(ratio)>10**-6:
        linEcc = linEcc - ratio
        func = linEcc-ecc*np.sin(linEcc)-Me
        deriv = 1-ecc*np.cos(linEcc)
        ratio = func/deriv

    trueAn = np.sqrt((1+ecc)/(1-ecc))*np.tan(linEcc/2)
    trueAn = 2*np.arctan(trueAn)

    return trueAn%(2*np.pi)

def trueToMean(tAn, ecc):
    """
    Return the mean anomaly based on the true anomaly

    INPUT:
        tAn : True anomaly float
        ecc : Eccentricity of the orbit float
    OUTPUT:
        The mean anomaly for the corresponding orbit
    """
    linEcc = 2*np.arctan(np.sqrt((1-ecc)/(1+ecc))*np.tan(tAn/2))
    return (linEcc-ecc*np.sin(linEcc))%(2*np.pi)

def S(x):
    """
    Stumpff S-function
    """
    if x < 0:
        y = (np.sinh(np.sqrt(-x))-np.sqrt(-x))/(np.sqrt(-x)**3)
    if x > 0:
        y = (np.sqrt(x)-np.sin(np.sqrt(x)))/(np.sqrt(x)**3)
    else:
        y = 1/6
    return y

def C(x):
    """
    Stumpff C-function
    """
    if x<0:
        y = (np.cosh(np.sqrt(-x))-1)/(-x)
    if x>0:
        y = (1-np.cos(np.sqrt(x)))/x
    else:
        y = 1/2
    return y

#TODO: Create fuzzy testing for lambert problem
def uLambert(r1,r2,dt0, primBody, nRev=0, DM=None):
    """
    Based on http://ccar.colorado.edu/imd/2015/documents/LambertHandout.pdf
    Solution to the Lambert problem
    The algorithm doesn't seem very robust and tends to fail for weird tof values

    INPUT:
        r1  : Initial position vector (x,y,z)
        r2  : Final position vector   (x,y,z)
        dt0 : Time of flight between r1 and r2 (in seconds)
        primBody : The primary body around which r1 and r2 are positioned, Body object
        DM : Optional, direction of motion, can take values of 1 and -1

    OUTPUT:
        r1  : Initial position vector (x,y,z)
        v1  : Initial velocity vector (x,y,z)
        r2  : Final position vector (x,y,z)
        v2  : Final velocity vector (x,y,z)
    """
    def _detCFunctions(psi):
        if psi > 1*10**-6:
            c2 = (1 - np.cos(np.sqrt(psi)))/psi
            c3 = (np.sqrt(psi) - np.sin(np.sqrt(psi)))/np.sqrt(psi**3)
        elif psi < -1*10**-6:
            c2 = (1 - np.cosh(np.sqrt(-psi)))/psi
            c3 = (np.sinh(np.sqrt(-psi))-np.sqrt(-psi))/np.sqrt((-psi)**3)
        else:
            c2 = 1/2
            c3 = 1/6
        return c2,c3

    def _detTOF(psi, A, mu, r1Mag, r2Mag):
        c2, c3 = _detCFunctions(psi)
        y = r1Mag+r2Mag+A*(psi*c3-1)/(np.sqrt(c2))
        chi = np.sqrt(y/c2)
        dt = (c3*np.power(chi,3)+A*np.sqrt(y))/np.sqrt(mu)
        return dt

    r1 = np.array(r1)
    r2 = np.array(r2)
    mu = primBody.mu

    if DM==None: # If we don't know the direction of motion we determine it
        nu1 = np.arctan2(r1[1], r1[0])
        nu2 = np.arctan2(r2[1], r2[0])
        dnu = (nu2 - nu1)%(2*np.pi)

        if 2*np.pi < dnu or dnu < 0:
            raise ValueError("dnu not in correct range")
        if dnu < np.pi:
            DM = 1
        else:
            DM = -1

    r1Mag = np.linalg.norm(r1)
    r2Mag = np.linalg.norm(r2)
    print(r1Mag)
    cosdnu = np.dot(r1,r2)/(r1Mag*r2Mag)

    A = DM*np.sqrt(r1Mag*r2Mag*(1+cosdnu))

    if A == 0:
        raise ValueError("Trajectory can't be computed")

    if nRev == 0:
        psi, psiUp, psiLow = 0, 4*np.pi**2, -4*np.pi
        c2, c3 = 1/2, 1/6
    else:
        psiUp, psiLow = 4*((nRev+1)**2)*np.pi**2, 4*(nRev**2)*np.pi**2
        psi = (psiUp+psiLow)/2
        c2, c3 = _detCFunctions(psi)

    #We verify that a solution exists for the time of flight and the nRev
    psiMini = fmin(_detTOF, psi, args=(A, mu, r1Mag, r2Mag))
    tofMini = _detTOF(psiMini[0], A, mu, r1Mag, r2Mag)
    if tofMini > dt0:
        raise ValueError("""
                        No solution : Time of flight too low

                        Number of revolutions   - {}
                        TOF Chosen              - {}
                        TOF Minimum for nRev    - {}
                        """.format(nRev, dt0, tofMini))
    dt = 0
    nIter = 0

    while abs(dt - dt0) > 1*10**-6:
        y = r1Mag+r2Mag+A*(psi*c3-1)/(np.sqrt(c2))
        if A > 0 and y < 0: #We make sure that y>0 by tweaking its value
            while y<0:
                c2,c3 = _detCFunctions(psi)
                psi += 0.1
                y = r1Mag+r2Mag+A*(psi*c3-1)/(np.sqrt(c2))


        chi = np.sqrt(y/c2)
        dt = (c3*np.power(chi,3)+A*np.sqrt(y))/np.sqrt(mu)

        if dt <= dt0:
            psiLow = psi
        else:
            psiUp = psi
        psi = (psiUp + psiLow)/2

        c2, c3 = _detCFunctions(psi)

        nIter+=1

        if nIter > 500:
            raise ValueError("""Maximum iterations reached in Lambert solver

            psi  - {}
            dt   - {}
            y    - {}
            dnu  - {}
            c2/3 - {} , {}

            Make sure your values are sensible and try again.
            """.format(psi, dt, y, dnu, c2, c3))

    f = 1 - y/r1Mag
    g = A*np.sqrt(y/mu)
    gDot = 1 - y/r2Mag
    v1 = (r2 - f*r1)/g
    v2 = (gDot*r2 - r1)/g

    return r1,v1,r2,v2

# def multiRevLambert(r1, r2, dt0, primBody, nRev=0, DM=None):


if __name__ == "__main__":
    from examples import Earth, Sun

    print(uLambert(np.array([5000,10000,2100]),np.array([-14600,2500,7000]),3600,primBody=Earth))
    print(uLambert(np.array([-26503064.373261325, 144693278.60351047, 119.32162712262743]),np.array([208034201.40482172, -1959743.5511186407, -5158244.7539406745]), 300*3600*24, primBody=Sun))
