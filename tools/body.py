class Body:
    def __init__(self, name=None, mu=0, r=0, primBody = None, Orbit=None, Atmosphere=None, Moons=None):
        self._name = name
        self._mu = mu
        self._r = r
        self._Orbit = Orbit
        self._Atmosphere = Atmosphere  #NOT_IMPLEMENTED
        self._Moons = Moons #NOT_IMPLEMENTED
        self._primBody = primBody

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if type(value) != str:
            raise TypeError("The body name must be a string")

        self._name = value

    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, value):
        if type(value) != int:
            raise TypeError("The gravitational parameter must be an integer")
        if value <= 0:
            raise ValueError("The gravitational parameter must be positive")
        self._mu = value

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        if type(value) != int:
            raise TypeError("The radius must be an integer")
        if value <= 0:
            raise ValueError("The radius must be positive")
        self._r = value

    @property
    def Orbit(self):
        return self._Orbit

    @Orbit.setter
    def Orbit(self, value):
        self._Orbit = value

    @property
    def primBody(self):
        return self._primBody

    @primBody.setter
    def primBody(self, value):
        if type(value) != Body:
            raise TypeError("The primary body must be another body")

        self._primBody = value
