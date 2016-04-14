"""
Module to do crude SN calculations for observations of a Lunar Occultation.

At the moment this is hardwired for La Palma/ULTRACAM/HiperCam
"""
import numpy as np
from astropy import units as u

# typical lunar sky background brightnesses for SDSS filters
# mag/sq arcsec
skyarr = np.array([ 13.6, 12.9, 12.6, 12.11, 11.83])
# extinction based on measured la palma values
extarr = np.array([0.50,0.19,0.09,0.05,0.04])
# zeropoints, based on WHT/ULTRACAM
# these zeropoints give 1 count per second
zparr = np.array([25.05, 26.88, 26.3, 26.33, 25.28])

gain = 1.2 #e-/ADU
rnoarr = [2,5] # guesses for rno in e-/pix for slow/fast

# how to scale counts for gtc from WHT. 
# scaled by ratio of apertures plus extra reflection
gtc_scalefactor = (730000./124700.) * 0.9

ps = 0.15 # arcsec/pix (guess based on VLT)

def filtSel(filt):
    name = filt.name
    filt_mapping = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'ha':2}
    return filt_mapping[name]
    
def getParams(filt, speed = 0):
    """Get parameters needed to calculate SN
    
    Parameters
    -----------
    filt : ~occult.filters.Filter
        filter of observations
    speed : int
        0 for slow, 1 for fast (default = 0)
        
    Returns
    -------
    sky, extinction, zp, rno
    """
    indx = filtSel(filt)
    zp = zparr[indx]
    # if it's Halpha, scale r-band zp appropriately
    # r has 27.4x the bandwidth
    # star must be 3.6 mags brighter
    if filt.name == 'ha':
        zp -= 3.6
    return skyarr[indx], extarr[indx], zp, rnoarr[speed]
        
@u.quantity_input(expT=u.s)
def sn(objmag, filt, expT, seeing, airmass, speed=0):
    expT = expT.to(u.s).value
    skymag, ext, zp, rno = getParams(filt, speed)
    objcounts = (10**((zp-objmag)/2.5)) * (10.**(-(ext*airmass)/2.5)) * expT
    objcounts *= gtc_scalefactor
    # number of pixels in seing disc
    pixels = np.pi * (seeing/ps)**2.
    
    skycounts = ((10**((zp-skymag)/2.5)) * (10.**(-(ext*airmass)/2.5)) *
                  expT * ps * ps)
    skycounts *= gtc_scalefactor
    
    sn = (objcounts/gain) / np.sqrt(objcounts*gain + 
                                    pixels*(skycounts*gain + rno**2))
    return sn  
                 
