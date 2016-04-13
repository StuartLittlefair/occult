from __future__ import (print_function)
import numpy as np
from scipy import special
from astropy import units as u, constants as const

from ..filters import sdss_filters

__all__ = ['make_fringe_pattern']

@u.quantity_input(vmoon=u.km/u.s)
def make_fringe_pattern(vmoon):
    """Compute the diffraction pattern from Lunar occultation by a point source
    
    The normalised diffraction pattern of a point source from a Lunar occulation
    is given by 
    
    .. math::
        
        I = I_0 \int_{-\infty}^x e^{\frac{i\piu^2}{D\lambda}} du
        
    where :math:`x` is the distance from the Moon's shadow limb, :math:`D` is the
    Lunar distance and :math:`\lambda` is the wavelength of observation.
        
    Parameters
    -----------
    vmoon :  ~astropy.units.Quantity
        Speed of lunar shadow, in km/second or equivalent units.
        
    Returns
    --------
    time : ~astropy.units.Quantity
        Array of times from start of occultation
    flux : np.ndarray
        Normalised flux from star
    """
    
    # moon distance
    # TODO: use time of occultation and jplehem to calculate actual moon distance
    moon_distance = 3.8e8 * u.m
    
    # number of points to calculate
    # TODO: improve to calculate at give time resolution with subdivision of bins
    nelem = 2400
    
    # number of wavelngths to calculate
    nwav = 200
    
    # Filter qualities (SDSS g')
    # TODO: pass in a filter object
    filt = sdss_filters['g']
    filter_start = filt.pivot - filt.fwhm/2
    filter_end   = filt.pivot + filt.fwhm/2
    wavelengths = u.AA*np.linspace(filter_start, filter_end, nwav)
    
    
    # define moon shadow positions from -40 to +60 meters
    x = u.m * np.linspace(-40, 60, nelem)
    time = (-x/vmoon).to(u.ms)
    
    # integral limits for Fresnel integral
    X, WAV = np.meshgrid(x,wavelengths)
    WLIMS = X * np.sqrt(np.pi/moon_distance/WAV)
    # WLIMS is dimensionless; just get values
    WLIMS = WLIMS.decompose().value
    
    """ 
    modfresnelp calulates
    
    .. math::
    
        \int_x^{\infty} e^{it^2} dt
     
    and returns a complex answer.   

    The calculation below transforms that into the integral we want, described
    in the docstring
    """
    fp, _ = special.modfresnelp(-WLIMS)
    F = np.absolute(-fp)
    
    # each row of F contains the diffraction pattern for a single wavelength.
    # Add these together and normalise
    F = F.sum(axis=0)
    F /= F.max()
    
    return time, F
 
       