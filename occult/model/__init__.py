from __future__ import (print_function)
import numpy as np
from scipy import special
from astropy import units as u, constants as const
from astropy.convolution import convolve

__all__ = ['make_fringe_pattern']

@u.quantity_input(moon_distance=u.km)
@u.quantity_input(vmoon=u.km/u.s)
@u.quantity_input(frame_rate=u.ms)
def make_fringe_pattern(moon_distance, vmoon, filt, source=None, telescope=None,
                        nelem=2400):
    """Compute the diffraction pattern from Lunar occultation by a given source.
    
    The normalised diffraction pattern of a point source from a Lunar occulation
    is given by 
    
    .. math::
        
        I = I_0 \int_{-\infty}^x e^{\frac{i\piu^2}{D\lambda}} du
        
    where :math:`x` is the distance from the Moon's shadow limb, :math:`D` is the
    Lunar distance and :math:`\lambda` is the wavelength of observation.
    
    For finite sources and telescopes, this point source diffraction pattern needs
    to be convolved with a 1D kernel. This 1D Kernel should be the chord the lunar
    limb makes through telescope aperture, or the projection of the sources intensity
    distribution from the lunar limb to the Earth 
    (see http://spiff.rit.edu/richmond/occult/bessel/bessel.html for more details).

        
    Parameters
    -----------
    moon_distance: ~astropy.units.Quantity
        Distance of Moon at time of occultation
    vmoon :  ~astropy.units.Quantity
        Speed of lunar shadow, in km/second or equivalent units.
    filt: ~occult.filter
        A Filter object representing the filter of observation
    source : ~occult.model.aperture.Aperture
        An Aperture object representing the occulted source
    telescope : ~occult.model.aperture.Aperture
        An Source object representing the occulted source
    nelem : integer
        Number of points to calculate, default=2400
    Returns
    --------
    time : ~astropy.units.Quantity
        Array of times from start of occultation
    flux : np.ndarray
        Normalised flux from star
    """
    # number of wavelngths to calculate
    nwav = 200
    
    # Filter qualities (SDSS g')
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
    
    # Now convolve with 1D Kernel from Source
    if source is not None:
        source.project(moon_distance)
        kernel = source.evaluate_1d_kernel(x)
        F = convolve(F, kernel, boundary='extend')
    
    # and telescope aperture
    if telescope is not None:
        kernel = telescope.evaluate_1d_kernel(x)
        F = convolve(F, kernel, boundary='extend')
        
    #normalise
    F /= F[-100:].mean()   
    return time, F
 

       
