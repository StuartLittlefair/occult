import numpy as np
import argparse

import occult.model as m
import occult.filters as f
import occult.model.apertures as ap
from occult.observing.rebin import rebin
from occult.observing.sn import sn
from astropy import units as u
from matplotlib import pyplot as plt
from scipy import stats

parser = argparse.ArgumentParser( description='Estimate resolution of Lunar Occultation')

parser.add_argument('--tel_diam', type=float, default = 10.2,
                    help='diameter of telescope in m : default 10.2')
parser.add_argument('--expT', type=float, default = 0.5,
                    help='exposure Time of observations in milliseconds : default=0.5')
parser.add_argument('--vmoon', type=float, default = 1,
                    help='speed of lunar shadow, in km/s : default 1')
parser.add_argument('--moon_dist', type=float, default = 384400,
                    help='distance to moon, in km : default 384,000')                    
parser.add_argument('--filter', default='g',
                    help='observations filter (u,g,r,i,z,ha) : default g')
parser.add_argument('--mag', default=8, type=float,
                    help='magnitude of star to observe : default 8')
parser.add_argument('--start', default=0.1, type=float,
                    help='initial resolution (mas) : default 0.1')
parser.add_argument('--end', default=2.0, type=float,
                    help='final resolution (mas) : default 2.0')
parser.add_argument('--step', default=0.2, type=float,
                    help='resolution step (mas) : default 0.2')
                                                           
args = parser.parse_args()

star = ap.UniformStar(1*u.R_sun, 1*u.pc)
tel = ap.Telescope(args.tel_diam*u.m)

moon_distance = args.moon_dist*u.km
vmoon = args.vmoon * u.km/u.s

if args.filter.lower() == 'ha':
    filt = f.ha
else:
    filt = f.sdss_filters[args.filter.lower()]

expT = args.expT*u.ms
# assume seeing 0.8 arcsecs, airmass 1.3, slow speed
signal_to_noise = sn(args.mag, filt, expT, 0.8, 1.3, 0)
print("Obtaining SNR of ", signal_to_noise)

# create a point source dataset at very high resolution
tfake, ffake = m.make_fringe_pattern(moon_distance, vmoon, filt, nelem=2800, telescope=tel)

# bin to correct time resolution
tmin = tfake.to(u.ms).value.min()
tmax = tfake.to(u.ms).value.max()
tbins = np.arange(tmin, tmax, args.expT)[1:-1]
tbin, fmodel = rebin(tbins, tfake.value, ffake)

resolutions = u.mas*np.arange(args.start, args.end, args.step)
nmonte = 50
dof = len(tbin)
for resolution in resolutions:
    
    star.angular_size = resolution
    num_rejected = 0
    for i in range(nmonte):
        tfake, ffake = m.make_fringe_pattern(moon_distance, vmoon, filt, nelem=2800, 
                                             source=star, telescope=tel)
        tbin, fdata = rebin(tbins, tfake.value, ffake)
        errs = fdata/signal_to_noise
        
        # add random noise
        fdata += np.random.normal(scale=errs)
        
        chisq = np.sum((fdata-fmodel)**2 / errs**2)
        pchisq = stats.chi2.sf(chisq, dof)
        # if < than 5% chance of dataset exceeding chisq, reject it
        if pchisq < 0.05:
            num_rejected += 1
    # have we rejected 95% of tests, if so we can resolve a source this size
    print(resolution, ' -> ', num_rejected / nmonte)
    if num_rejected / nmonte > 0.95:
        break
        
print('Estimated resolution is %.2f mas' % resolution.value)
chisq = np.sum((fdata-fmodel)**2 / errs**2)
pchisq = stats.chi2.sf(chisq, dof)
print('Chisq = %f with %d DOF (P=%f percent)' % (chisq,dof,100*pchisq))
plt.style.use('ggplot')
plt.plot(tbin,fmodel,label='Point Source')
plt.errorbar(tbin,fdata,yerr=errs,fmt='.')
plt.title('Estimated resolution is %.2f mas\nChisq = %f with %d DOF (P=%f percent)' 
           % (resolution.value,chisq,dof,100*pchisq))
plt.show()