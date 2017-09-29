from orbits.orbit import Orbit
from examples import Earth, Mars, Sun
from plot.plotlyManager import PlotManager3D
import numpy as np
from tools.trajOpt import getC3

if __name__ == '__main__':
    initialShift = 500
    Earth.Orbit.updTime(initialShift*3600*24)
    Mars.Orbit.updTime(initialShift*3600*24)
    pointDict = {   'color' : 'rgba(0,0,0,255)',
                    'size'  : 5}
    orbitDict = {   'color' : 'rgba(100, 100, 100, 255)',
                    'width' : 2}
    pltMgr = PlotManager3D(iterOrbit = 300,
    iterAnim = 100, title = "Earth -> Mars in 168 days",
    range = Mars.Orbit.getApoapsis()*1.5, pointBase = pointDict,
    orbitBase = orbitDict, filename="MarsTransfer.html")

    days = 168*5

    pltMgr.addBody(Sun)
    pltMgr.addOrbit(Earth.Orbit, paramDict = {'color' : 'blue'})
    pltMgr.addOrbit(Mars.Orbit, paramDict = {'color' : 'orange'})
    pltMgr.animPoint(Earth.Orbit, tAnStart=Earth.Orbit.tAn, time = 3600*24*days)
    pltMgr.animPoint(Mars.Orbit, tAnStart=Mars.Orbit.tAn, time = 3600*24*days)


    Mars.Orbit.updTime(3600*24*days)
    transferOrbit = Orbit.fromLambert(Earth.Orbit.cartesianCoordinates(),
                Mars.Orbit.cartesianCoordinates(), 3600*24*days, primBody=Sun, nRev=2,
                    name="Transfer Orbit {} days".format(days))


    tTrans1 = transferOrbit.tAn
    tTrans2 = transferOrbit.tAnAtTime(3600*24*days)

    pltMgr.addOrbit(transferOrbit, paramDict={'color' : 'green'})
    pltMgr.animPoint(transferOrbit, tAnStart=transferOrbit.tAn, time = 3600*24*days)
    # print(getC3(Earth, Mars.Orbit.cartesianCoordinates(), days*3600*24))
    pltMgr.show()
