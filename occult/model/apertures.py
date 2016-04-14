from abc import ABCMeta, abstractmethod

import numpy as np
from astropy import units as u

class Aperture:
    __metaclass__ = ABCMeta
    @abstractmethod
    def evaluate_1d_kernel(self, x): pass
    
class Disk(Aperture):
    """
    An aperture representing a uniform, circular Disk.
    
    Such an aperture can be subclassed for e.g. a simple star, or a telescope aperture
    
    Every aperture can return a 1D Kernel used for convolving the point source
    occultation profile.
    """
    @u.quantity_input(radius=u.m)
    def __init__(self,radius):
        self.radius = radius.to(u.m)
        # scale is used when projecting onto the ground
        self.scale  = 1.0
    
    def evaluate_1d_kernel(self, x):
        """
        Return array of 1d kernel, evaluated on the same grid as x.
        """
        R = (self.radius*self.scale)
        inner = R**2 - np.fabs(R-x)**2
        mask = inner >= 0
        kernel = np.sqrt(inner[mask]).value
        if len(kernel)%2 == 0:
             #kernel must have odd dimensions
             kernel = np.append(kernel,0)
        return kernel

class Telescope(Disk):
    @u.quantity_input(diameter=u.m)
    def __init__(self, diameter):
        super(Telescope, self).__init__(diameter/2)
    
class UniformStar(Disk):
    """
    An aperture for modelling the effects of a finite, uniform stellar disk
    """
    @u.quantity_input(radius=u.R_sun)
    @u.quantity_input(distance=u.pc)
    def __init__(self, radius, distance):
        super(UniformStar, self).__init__(radius)
        self.distance = distance
        
    @property
    def angular_size(self):
        return ((self.radius/self.distance).decompose() * u.rad).to(u.mas)
        
    @angular_size.setter
    @u.quantity_input(angular_size=u.rad)
    def angular_size(self, angular_size):
        # consider radius as fixed, move distance to get right angular size
        self.distance = self.radius/angular_size.to(u.rad).value

    @u.quantity_input(moon_distance=u.m)
    def project(self, moon_distance):
        """Calculate appropriate scale for a lunar occultation"""
        projected_size = moon_distance * np.tan(self.angular_size)
        self.scale = projected_size / self.radius