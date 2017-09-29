import tools.utils
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ITERORBIT = 200

class PlotManager():
    def __init__(self):
        self.orbitList = {}
        self.fig, self.ax = plt.subplots()
        plt.axis('equal')

    def addOrbit(self, Orbit):
        self.orbitList[Orbit.id] = Orbit

    def plotOrbit(self, orbitID):
        Orbit2D = self.orbitList[orbitID]

        tAnArray = np.linspace(0, 2*np.pi, ITERORBIT)
        xArray, yArray = np.array([]), np.array([])

        for i in np.nditer(tAnArray):
            # i = utils.meanToTrue(i, Orbit2D.e)
            x,y,_ = Orbit2D.cartesianCoordinates(i)

            xArray = np.append(xArray,x)
            yArray = np.append(yArray,y)

        self.ax.plot(xArray, yArray)


        # ELLIPSE METHOD, far more efficient, need to debug however
        #BUG: Doesnt'display properly based on aPe, fix the ellipse patches


        # SMA = Orbit2D.a #SMajor
        # SMI = Orbit2D.a*np.sqrt(1-Orbit2D.e**2)
        # linEccentricity = np.sqrt(SMA**2-SMI**2)
        # xF, yF = utils.sphereToCartesian2D(-Orbit2D.aPe, linEccentricity)
        # body = mpatches.Ellipse([xF, yF], SMA*2, SMi*2, Orbit2D.aPe, fill=False, edgecolor=np.random.rand(3), linestyle='--')
        # self.ax.add_patch(body)

    def plotPoint(self, orbitID, area=25):
        Orbit2D = self.orbitList[orbitID]
        x, y, _ = Orbit2D.cartesianCoordinates(Orbit2D.tAn)

        self.ax.scatter(x,y, s=area)

    def plotBody(self, body):
        xArray = np.array([])
        yArray = np.array([])

        radius = body.r

        body = mpatches.Circle([0,0],radius)
        self.ax.add_patch(body)

    def show(self):
        self.ax.set_aspect('equal')
        self.ax.autoscale()
        plt.show()
