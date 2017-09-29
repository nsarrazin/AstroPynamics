import tools.utils
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d


ITERORBIT = 200

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().

      Thanks to karlo for this function : https://stackoverflow.com/a/31364297/7886572
    '''
    print(ax.get_xlim3d())
    print(ax.get_ylim3d())
    print(ax.get_zlim3d())
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


class PlotManager3D():
    def __init__(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')

        self.ax.set_xlabel('X axis')
        self.ax.set_ylabel('Y axis')
        self.ax.set_zlabel('Z axis')
        self.orbitList = {}

    def addOrbit(self, Orbit):
        self.orbitList[Orbit.id] = Orbit

    def plotBody(self, body):
        r = body.r

        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = r * np.outer(np.cos(u), np.sin(v))
        y = r * np.outer(np.sin(u), np.sin(v))
        z = r * np.outer(np.ones(np.size(u)), np.cos(v))
        print(np.amax(x))
        print(np.amax(y))
        print(np.amax(z))
        self.ax.plot_surface(x, y, z, color='b')

    def plotOrbit(self, orbitID, color=None):
        Orbit3D = self.orbitList[orbitID]
        xArray, yArray, zArray = np.array([]), np.array([]), np.array([])
        tAnArray = np.linspace(0, 2*np.pi, ITERORBIT)

        for i in np.nditer(tAnArray):
            i = tools.utils.meanToTrue(i, Orbit3D.e)
            x,y,z = Orbit3D.cartesianCoordinates(i)
            xArray = np.append(xArray, x)
            yArray = np.append(yArray, y)
            zArray = np.append(zArray, z)

        self.ax.plot(xArray, yArray, zArray, c=color, lw = 0.5)

    def show(self):
        # self.ax.set_aspect('equal')
        # self.ax.axis('equal')
        # set_axes_equal(self.ax)
        # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        # x = np.cos(u)*np.sin(v)
        # y = np.sin(u)*np.sin(v)
        # z = np.cos(v)
        # self.ax.plot_wireframe(x, y, z, color="r")
        self.ax.view_init(elev = 90, azim= 0)
        plt.show()
