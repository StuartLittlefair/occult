import numpy as np
def rebin(xbins,x,y):
    digitized = np.digitize(x,xbins)
    xbin = []
    ybin = []
    for i in range(0,len(xbins)):
            bin_y_vals = y[digitized == i]
            bin_x_vals = x[digitized == i]
            xbin.append(bin_x_vals.mean())
            ybin.append(bin_y_vals.mean())
    xbin = np.array(xbin)
    ybin = np.array(ybin)
    return (xbin,ybin)