from astropy import units as u
footcandle = u.lumen/u.imperial.foot**2
lambert = 1*u.cd/u.cm/u.cm
nanolambert = 1.0e-9*lambert

@u.quantity_input(cusp_angle=u.deg)
def lunar_background(sunlit_fraction, airmass, cusp_angle, ext=0.172):
    """
    Estimate lunar V-band background, based on Schaefer, Bulder & Bourgeois 1992.
    
    Parameters
    ----------
    sunlit_fraction : float
        Fraction of Lunar disc which is illuminated (1=Full Moon)
    airmass : float
        Airmass of observations
    cusp_angle : ~astropy.units.Quantity
        Location of star on Lunar limb during occultation
        
    Returns
    --------
    ugriz : ~astropy.units.Quantity
        ugriz background brightnesses in Mags/arcsec**2
    """
    # Illuminance of Sun
    # all illuminances in footcandles!
    Isun = 1.174e4 

    # covert sunlit fraction to phase angle
    alpha = np.arccos(2*sunlit_fraction-1) * u.rad
    alpha = alpha.to(u.deg).value
    
    # visual magnitude of entire moon.
    m = -12.73 + 0.026*np.fabs(alpha) + 4.0e-9*alpha**4
    # Illuminance of Moon
    Imoon = 10**(-0.4*(m+16.57)) 
    
    # Calculate the complement for calculating Earthshine
    alpha = 180-np.fabs(alpha)
    m = -12.73 + 0.026*np.fabs(alpha) + 4.0e-9*alpha**4
    Imoon_prime = 10**(-0.4*(m+16.57)) 

    # Illuminance of Earth as viewed from the Moon, found
    # by scaling Moon's brightness by area and Albedos
    Iearth = 78*Imoon_prime
    
    # Apparent surface brightness of Earthshine
    # all surface brightnesses in nanolamberts!
    Bmoon0 = 1.65e9 * 10.0**(-0.4*ext*airmass)
    B_es = Bmoon0*Iearth/Isun
        
    # Now the scattered moonlight
    theta_moon = 0.26 # angular size of moon, degrees
    # effective separation between source and center of Moon
    theta = theta_moon * (1.0 - 0.4*np.exp(-cusp_angle.to(u.deg).value/30))
    theta = theta * (np.cos(cusp_angle)**2 + (1-sunlit_fraction+np.sin(cusp_angle))**2)**0.5

    B_sc = 6.25e7 * Imoon / theta**2.0 
    B_sc = B_sc * (10.0**(-0.4*ext*airmass) - 10.0**(-0.8*ext*airmass)) 
    B_total = B_sc+B_es

    # now convert from nanolamberts to mag/sq arcsec
    V = (20.7233-np.log(B_total.value/34.08))/0.92104
    U = V+0.2
    B = V+0.5
    R = V-0.3
    I = V-1.0
    magsq = np.array([U,B,V,R,I])
    return u.Magnitude(ubv2ugr(magsq))/u.arcsec**2
    
def ubv2ugr(ubv):
    """convert UBVRI to ugriz using the Smith et al equations"""
    u,b,v,r,i = ubv
    g_sdss = v+(0.54*(b-v))-0.07
    r_sdss = v-(0.44*(b-v))+0.12
    u_sdss = (1.33*(u-b))+1.12+g_sdss
    i_sdss = r_sdss-(1.00*(r-i))+0.21
    z_sdss = r_sdss-(1.65*(r-i))+0.38
    return np.array([u_sdss, g_sdss, r_sdss, i_sdss, z_sdss])