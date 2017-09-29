import numpy as np
from astropy.time import Time, TimeDelta
from tools import utils

class OrbitalElements:
    def __init__(self, SMA=None, eccentricity=0, inclination=0, lAn=0, tAn=0, primBody=None, aPe = 0, epoch=Time('2000-01-01'),
                name=None):
        self._a = SMA #[0-inf] km
        self._e = eccentricity #[0-1]
        self._i = inclination
        self._lAn = lAn
        self._aPe = aPe #[0-359] deg
        self._tAn = tAn #000
        self._primBody = primBody
        self._id = name
        self._epoch = epoch

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self,SMA):
        if SMA<=0:
            raise ValueError("Semi-Major Axis inferior or equal to zero")
        self._a = SMA

    @property
    def e(self):
        return self._e

    @e.setter
    def e(self,eccentricity):
        if 1<=eccentricity<0:
            raise ValueError("eccentricity not between 0 and 1")
        self._e = eccentricity

    @property
    def i(self):
        return self._i

    @i.setter
    def i(self, value):
        if value<0:
            raise ValueError("negative inclination")
        self._i = value%2*np.pi

    @property
    def lAn(self):
        return self._lAn

    @lAn.setter
    def lAn(self, value):
        if value<0:
            raise ValueError("negative longitude of the ascending node")
        self._lAn = value%2*np.pi

    @property
    def aPe(self):
        return self._aPe

    @aPe.setter
    def aPe(self, value):
        if value<0:
            raise ValueError("negative argument of perapsis")
        self._aPe = value%2*np.pi

    @property
    def tAn(self):
        return self._tAn

    @tAn.setter
    def tAn(self, value):
        if value<0:
            raise ValueError("negative true anomaly")
        self._tAn = value%(2*np.pi)

    @property
    def epoch(self):
        return self._epoch

    @epoch.setter
    def epoch(self, value):
        if not isinstance(value, Time):
            raise TypeError("The epoch should be an astropy Time object")
        self._epoch = value

    @property
    def primBody(self):
        return self._primBody

    @primBody.setter
    def primBody(self,value):
        if type(value) != dict:
            raise TypeError("The body should be a dict, check bodies.py for ref")
        self._primBody = value

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        if type(value) != str:
            raise TypeError("Name provided is not a string")
        self._id = value

class Orbit(OrbitalElements):
    def __init__(self,  SMA = None, eccentricity=0, inclination=0, lAn=0, tAn=0, primBody=None, aPe = 0, epoch=Time('2000-01-01'),
                name=None):
        OrbitalElements.__init__(self,  SMA, eccentricity, inclination, lAn, tAn, primBody, aPe, epoch,
                    name)
    @classmethod
    def fromElements(cls, semiMajorAxis=None, eccentricity = 0, inclination=0,
                    lAn=0, primBody = None, argumentPeriapsis = 0,
                    trueAnomaly = 0, refEpoch=Time('2000-01-01'), name=None):
        """
        fromElements is callable from the class itself in order to create an
        orbit from the classical Keplerian elements.
        """

        if semiMajorAxis==None:
            raise ValueError("Semi-Major Axis required !")

        if name == None:
            raise ValueError("The orbit needs a name !")

        return cls(SMA=semiMajorAxis, eccentricity=eccentricity, inclination=inclination, primBody=primBody, aPe=argumentPeriapsis, tAn=trueAnomaly, lAn=lAn, epoch=refEpoch, name=name)

    @classmethod
    def fromApsis(cls, periapsis, apoapsis, inclination=0, lAn=0,
                    primBody = None, argumentPeriapsis = 0, trueAnomaly = 0, refEpoch=Time('2000-01-01'),
                    name=None):
        """
        fromApsis is to be called from the class itself as well, replaces the
        need for SMA and eccentricity by Apoapsis and Periapss
        """
        if name == None:
            raise ValueError("The orbit needs a name !")

        SMA = (periapsis+apoapsis)/2
        eccentricity = (apoapsis-periapsis)/(apoapsis+periapsis)

        return cls(SMA=SMA, eccentricity=eccentricity, inclination=inclination, primBody=primBody, aPe=argumentPeriapsis, tAn=trueAnomaly, lAn=lAn, epoch=refEpoch, name=name)

    @classmethod
    def fromStateVector(cls, r1, v1, primBody = None, refEpoch=Time('2000-01-01'),
                    name=None):
        """
        fromStateVector determines orbital parameters based on position & velocity vectors

        the reference epoch is associated with the state vector
        """
        if name == None:
            raise ValueError("The orbit needs a name !")

        mu = primBody.mu

        r = np.linalg.norm(r1)
        v = np.linalg.norm(v1)

        vr = np.dot(r1,v1)/r
        h = np.cross(r1,v1)
        hMag = np.linalg.norm(h)

        i = np.arccos(h[2]/hMag)

        N = np.cross([0,0,1], h)
        NMag = np.linalg.norm(N)

        lAn = np.arccos(N[0]/NMag)
        if N[1]<0:
            lAn = 2*np.pi - lAn

        e = 1/mu*((v**2-mu/r)*r1-r*vr*v1)
        eMag = np.linalg.norm(e)

        if 0 > eMag >= 1:
            raise ValueError("e = "+str(e)+", not an ellipse")

        aPe = np.arccos(np.dot(N/NMag, e/eMag))
        if e[2] < 0 :
            aPe = 2*np.pi - aPe

        tAn = np.arccos(1/eMag*((hMag**2/(mu*r))-1))
        if vr < 0:
            tAn = 2*np.pi-tAn

        SMA = -(mu/(2*(v**2/2 - mu/r)))
        return cls(SMA=SMA, eccentricity=eMag, inclination=i, primBody=primBody, aPe=aPe, tAn=tAn, lAn=lAn, epoch=refEpoch, name=name)

    @classmethod
    def fromLambert(cls,r1,r2,t, DM=None, nRev=0, primBody = None, refEpoch=Time('2000-01-01'),name=None):
        """
        fromLambert solves the lambert problem and use the first state vector
        returned to determine orbital parameters
        """
        try:
            r1,v1,_,_ = utils.uLambert(r1,r2,t,primBody,nRev=nRev, DM=DM)
        except RuntimeWarning:
            print("ERROR : Runtime warning detected, incorrect values ?")
            return None
        return Orbit.fromStateVector(r1, v1, primBody, refEpoch, name)

    def getRadius(self,trueAnomaly):
        return (self.a*(1-self.e**2))/(1+self.e*np.cos(trueAnomaly))

    def currentAltitude(self):
        return self.getRadius(self.tAn)

    def altitudeFromRadius(self,radius):
        alt = radius-self.primBody.r
        if alt < 0:
            raise ValueError("Altitude smaller than zero !")
        return alt

    def getApoapsis(self):
        apAngle = (self.aPe + 180)%(2*np.pi)
        return self.getRadius(apAngle)

    def getPeriapsis(self):
        return self.getRadius(self.aPe)

    def getPeriod(self):
        return (2*np.pi)*np.sqrt(self.a**3/self.primBody.mu)

    def getMeanMotion(self):
        return ((2*np.pi)/(self.getPeriod()))

    def tAnAtTime(self, sec):
        actualMe = utils.trueToMean(self.tAn, self.e)
        updMe = actualMe+self.getMeanMotion()*sec
        updtAn = utils.meanToTrue(updMe, self.e)
        return updtAn

    def updTime(self, sec):
        updtAn = self.tAnAtTime(sec)
        self.tAn = updtAn

    def isCrashing(self):
        try:
            altPe = self.altitudeFromRadius(self.getPeriapsis())
        except ValueError:
            return True
        return False

    def cartesianCoordinates(self, trueAnomaly=None):
        lAn = self.lAn
        aPe = self.aPe

        if trueAnomaly == None:
            tAn = self.tAn
        else:
            tAn = trueAnomaly

        i = self.i
        r = self.getRadius(tAn)

        x = r * ( np.cos(lAn)*np.cos(aPe+tAn) - np.sin(lAn)*np.sin(aPe+tAn)*np.cos(i))
        y = r * ( np.sin(lAn)*np.cos(aPe+tAn) + np.cos(lAn)*np.sin(aPe+tAn)*np.cos(i))
        z = r * ( np.sin(i)*np.sin(aPe+tAn))

        return x,y,z

    def getSpeed(self, trueAnomaly=None):
        if trueAnomaly == None:
            tAn = self.tAn
        else:
            tAn = trueAnomaly

        r = self.getRadius(tAn)
        mu = self.primBody.mu

        return np.sqrt(mu*(2/r -1/self.a))
