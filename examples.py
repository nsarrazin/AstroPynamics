from tools.body import Body
from orbits.orbit import Orbit
import numpy as np

Sun = Body(name="Sun",mu = 1.327124400*10**11, r = 695700)

_earthOrbit = Orbit.fromElements(semiMajorAxis=149597887.1558,
eccentricity=0.01671022, inclination = np.radians(0.00005),
lAn = np.radians(348.7394), primBody=Sun,
argumentPeriapsis=np.radians(114.2078),
trueAnomaly=6.23837308813, name="Earth Orbit")

Earth = Body(name="Earth", mu = 398600.4418, r = 6371,
Orbit=_earthOrbit, primBody=Sun)

_marsOrbit = Orbit.fromElements(semiMajorAxis=227936637.2418,
eccentricity=0.09341233, inclination = np.radians(1.85061),
lAn = np.radians(49.57854), primBody=Sun,
argumentPeriapsis=np.radians(286.4623),
trueAnomaly=0.40848952017, name="Mars Orbit")

Mars = Body(name="Mars", mu = 4.2828372854187757e+04, r = 3390, Orbit=_marsOrbit, primBody=Sun)
