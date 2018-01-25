import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sb
import spectroseti.apf as apf
import spectroseti.utilities as util

import scipy.signal as sg


from scipy.interpolate import CubicSpline
from sklearn.neighbors import KernelDensity
from sklearn.cluster import MeanShift, estimate_bandwidth

def meanshift_trend(spec,order):
    bstar1 = fits.open('/media/nate/DATA/Spectra/apf_bstar/ralh.272.fits')[0].data
    # shoudl take an order, deblaze it(B star), then do a meanshift.
    # then it should turn the main cluster into a continuum, meadian fit this, divide by it
    corrector_value = 100
    bstarfixed = (spec[order, :]+1) / util.savitzky_golay(sg.medfilt(bstar1[order, :]+1, kernel_size=501), 51, 4)

    bandwidth = estimate_bandwidth(bstarfixed[:,np.newaxis], quantile=0.1)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(bstarfixed[:,np.newaxis])
    xvals = np.arange(4608)
    test = np.array(bstarfixed)
    labels = ms.labels_
    test[labels != 0] = np.median(test[labels == 0])

    med_test = sg.medfilt(test, kernel_size=101)

    return bstarfixed/med_test

# I think this works!

# Incorporate it into the deblazig code, leave meanshift threasholding as is, and see if it imporves it.
# then we can use the B Stars for posterity