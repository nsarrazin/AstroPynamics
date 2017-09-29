import numpy as np
from orbits.orbit import Orbit
from examples import Earth, Sun
from plot.plotlyManager import PlotManager3D
from tools.utils import uLambert

nTest = 100


pointDict = {   'color' : 'rgba(0,0,0,255)',
                'size'  : 5}
orbitDict = {   'color' : 'rgba(100, 100, 100, 255)',
                'width' : 2}

pltMgr = PlotManager3D(iterOrbit = 300,
iterAnim = 100, title = "Lambert test",
range = 1E6, pointBase = pointDict,
orbitBase = orbitDict, filename="LambertTest.html")

for j in range(nTest):
    if j%100 == 0:
        print(j)

    hPa = np.random.uniform(low = 150, high = 100*10**4)
    hPe = np.random.uniform(low = 150, high = 100*10**4)
    i = np.random.uniform(low=0, high=np.pi)
    aPe = np.random.uniform(high = 2*np.pi)
    lAn = np.random.uniform(high = 2*np.pi)
    tAn = np.random.uniform(high = 2*np.pi)

    if hPe > hPa:
        hPe, hPa = hPa, hPe

    orbitTest = Orbit.fromApsis(hPe, hPa, inclination=i, argumentPeriapsis=aPe, lAn = lAn, trueAnomaly=tAn, name="test orbit - {}".format(j), primBody = Earth)

    pltMgr.addOrbit(orbitTest)
    # plt.plotOrbit(orbitTest.id)
    #

    r1 = orbitTest.cartesianCoordinates()

    tof = np.random.uniform()*orbitTest.getPeriod()*0.1
    # tof=3600*24
    orbitTest.updTime(tof)

    r2 = orbitTest.cartesianCoordinates()
    try:
        orbitComparison = Orbit.fromLambert(r1, r2, tof, primBody = Earth)
    except:
        print("""ERROR LAMBERT SOLVER
        KEPLERIAN ELEMENTS INITIAL ORBIT:
                hPa   - {}
                hPe   - {}
                i   - {}
                aPe - {}
                lAn - {}
                tAn - {}
                tof - {}""".format(hPa,hPe,i,aPe,lAn,tAn,tof))

        continue
    print("""COMPARISON OF KEPLERIAN ELEMENTS :
            a   - {}
            e   - {}
            i   - {}
            aPe - {}
            lAn - {}
            tAn - {}
            """.format(a-orbitComparison.a, e-orbitComparison.e, i-orbitComparison.i, aPe-orbitComparison.aPe, lAn-orbitComparison.lAn, tAn - orbitComparison.lAn))
pltMgr.show()

# testOrbit = Orbit.fromLambert(np.array([5000,10000,2100]),np.array([-14600,2500,7000]),3600,primBody=Earth, name="test")
# print(testOrbit.cartesianCoordinates())
# testOrbit.updTime(3600)
# print(testOrbit.cartesianCoordinates())
# print(uLambert(np.array([-26503064.373261325, 144693278.60351047, 119.32162712262743]),np.array([208034201.40482172, -1959743.5511186407, -5158244.7539406745]), 300*3600*24, primBody=Sun))
