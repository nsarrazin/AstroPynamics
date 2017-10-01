# AstroPynamics
A pure Python library to solve common orbital mechanics problem

## Orbit definition
### Basic class methods
The Orbit class includes a few class methods which allows for fast and easy Orbit generation.
Regardless of the input elements all orbits requires the following inputs :
* primBody  = a Body object used to determine the gravitational parameter
* name      = Every orbit created needs a name
#### fromElements
This uses the basic set of keplerian elements :
* Semi-major Axis (a)
* Eccentricity (e)
* Inclination (i)
* Longitude of the ascending node (lAn)
* Argument of Periapsis (aPe)
* True Anomaly (tAn)
```python
from orbits.orbit import Orbit
from examples import Sun
import numpy as np
earthOrbit = Orbit.fromElements(a=149597887.1558, e=0.01671022, i = np.radians(0.00005), lAn = np.radians(348.7394), primBody=Sun, aPe=np.radians(114.2078), tAn=6.23837308813, name="Earth Orbit")

print(earthOrbit)
```
Which prints :
```
    ID       - Earth Orbit
    primBody - Sun

    a    - 149597887.1558
    e    - 0.01671022
    i    - 8.726646259971648e-07
    lAn  - 6.086650761429513
    aPe  - 1.99330214145918
    tAn  - 6.23837308813

    epoch - 2000-01-01 00:00:00.000
````
#### fromApsis
fromApsis replaces the semimajor axis and the eccentricity by the height of the Apoapsis (hPa) and Periapsis (hPe)
```python
marsOrbit = Orbit.fromApsis(hPa = 247644270.465, hPe=220120745.788,  i = np.radians(1.85061), lAn = np.radians(49.57854), primBody=Sun, aPe=np.radians(286.4623), tAn=0.40848952017, name="Mars Orbit")

print(marsOrbit)
```
Which prints :
```
    ID       - Mars Orbit
    primBody - Sun

    a    - 233882508.1265
    e    - 0.05884049409567938
    i    - 0.03229923767033226
    lAn  - 0.8653087613317094
    aPe  - 4.999710317835753
    tAn  - 0.40848952017

    epoch - 2000-01-01 00:00:00.000
```
#### fromStateVector
This replaces the set of keplerian elements by a position 3D vector and a velocity 3D vector
```python
from examples import Earth

testOrbit = Orbit.fromStateVector((-6045, -3490, 2500), (-3.457, 6.618, 2.533), primBody=Earth, name = "State Vector Demo")

print(testOrbit)
```
Which prints :
```
    ID       - State Vector Demo
    primBody - Earth

    a    - 8788.081767279667
    e    - 0.17121118195416898
    i    - 2.6747036137846094
    lAn  - 4.455464041223287
    aPe  - 0.35025511728002945
    tAn  - 0.496472955354359

    epoch - 2000-01-01 00:00:00.000
```
### Universal multi-rev Lambert solver
Orbits can also be generated using a multi revolution universal Lambert solver. It uses the following inputs:

* r1 - 3D position vector
* r2 - 3D position vector
* tof - Time of flight between r1 and r2 in seconds

And the following are optional :

* nRev - The number of full revolution completed (default is 0)
* DM - The direction of motion, if left empty the solver will determine the optimal one

```python
lambertDemo = Orbit.fromLambert((5000,10000,2100), (-14600,2500,7000), 3600, primBody=Earth, name="Lambert Demo")
print(lambertDemo)
```
Which will output :
```
ID       - Lambert Demo
primBody - Earth

a    - 20002.884935993607
e    - 0.4334874513211504
i    - 0.5269331332631371
lAn  - 0.7784202841672524
aPe  - 0.535923312928038
tAn  - 6.123135458171056

epoch - 2000-01-01 00:00:00.000
```
The universal solver also supports hyperbolic trajectory :
```python
from orbits.hyperbola import Hyperbola
lambertDemo = Hyperbola.fromLambert((5000,10000,2100), (-14600,2500,7000), 1500, primBody=Earth, name="Lambert Hyperbola")
print("SMA - {}, eccentricity - {}".format(lambertDemo.a, lambertDemo.e))
```
Will output :
```
SMA - -2918.7342956582274, eccentricity - 4.215326772456562
```
It also supports multi-rev generation :
![0-rev](https://i.imgur.com/LWYBQ2X.gif)
![1-rev](https://i.imgur.com/ZXUqvlM.gif)
![2-rev](https://i.imgur.com/lpLI3nQ.gif)
## Plotting
### plotlyManager
TODO
### Plotting body
TODO
### Plotting orbits
TODO
### Plotting points
TODO
### Animation
TODO
## Maneuvers

## Trajectory optimization
### getC3

### Porkchop Generator
