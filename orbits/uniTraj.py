import numpy as np
from astropy.time import Time, TimeDelta
import utils

class TrajectoryElements:
    """
    TODO: Implement properties getter/setter and sanity checks
    """
    def __init__(self, universalAnomaly, radiusAtZero, velocityAtZero, reciprocalSMA, primaryBody, epoch, name):
        self._uAnom = universalAnomaly
        self._r0 = radiusAtZero
        self._v0 = velocityAtZero
        self._alpha = reciprocalSMA
        self._primBody = primaryBody
        self._epoch = epoch
        self._id = name

    @property
    def uAnom(self):
        return self._uAnom

    @uAnom.setter
    def uAnom(self, uAnomaly):
        self._uAnom = uAnomaly

    @property
    def r0(self):
        return self._r0

    @r0.setter
    def r0(self, radiusAtZero):
        if len(radiusAtZero) != 3 or type(radiusAtZero) != tuple:
                raise ValueError("radiusAtZero should be a tuple of length 3")
        self._r0 = radiusAtZero

    @property
    def v0(self):
        return self._v0

    @v0.setter
    def v0(self, velocityAtZero):
        if len(velocityAtZero) != 3 or type(velocityAtZero) != tuple:
                raise ValueError("velocityAtZero should be a tuple of length 3")
        self._v0 = velocityAtZero

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, reciprocalSMA):
        if type(reciprocalSMA) != int and type(reciprocalSMA) != float:
            raise TypeError("The reciprocal of the semi major axis should be an integer or a float")
        self._alpha = reciprocalSMA

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

class Trajectory(TrajectoryElements):
    # TODO: Implement universal kepler solver
    # TODO: Implement lagrange coefficients calculation
    # TODO: Allow conversion between Orbit class and Trajectory for e<1
    def __init__(self, universalAnomaly, radiusAtZero, velocityAtZero, reciprocalSMA, primaryBody, epoch, name):
        TrajectoryElements.__init__(self, universalAnomaly, radiusAtZero, velocityAtZero, reciprocalSMA, primaryBody, epoch, name)
    #
    # def uKeplerSolver(self, dt):
    #     nMax = 1000
    #     mu = self.primBody["gParam"]
    #     uAnom = np.sqrt(mu)*np.abs(self.alpha)*dt
    #     h = np.linalg.norm(np.cross(self.r0, self.v0))
    #     ratio = 1
    #     n = 0
    #     while ratio < 10**-8 and n <= nMax:
    #         C=utils.stumpC(self.alpha*uAnom**2)
    #         S=utils.stumpS(self.alpha*uAnom**2)
    #
    #         r0Norm = np.linalg.norm(self.r0)
    #
    #         f = (r0Norm*self.v0[0])/np.sqrt(mu)*uAnom**2*C + (1-self.alpha*r0Norm)*uAnom**3*S + r0Norm*uAnom-sqrt(mu)*dt
    #         df =
