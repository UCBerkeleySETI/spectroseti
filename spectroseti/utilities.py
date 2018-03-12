import math

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sg
from astropy.io import fits
from scipy.interpolate import CubicSpline
from scipy.ndimage.filters import percentile_filter
from scipy.signal import convolve2d
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.neighbors import KernelDensity

__author__ = 'nate'


def cantor(a, b):
    """Cantor pairing function, used to give unique int name to each observation"""
    a = int(a)
    b = int(b)
    return (a + b) * (a + b + 1) / 2 + b


def decantor(z):
    """Inverse Cantor"""
    w = math.floor(math.sqrt(8 * z + 1) / 2 - 0.5)
    t = ((w + 1) * w) / 2
    y = z - t
    x = w - y
    return int(x), int(y)


def filterbypercentile(arr, top, bottom):
    topp = np.percentile(arr, top)
    bottomp = np.percentile(arr, bottom)
    shortened = [bottomp < x < topp for x in arr]
    return shortened


def reject_outliers(data, m=2.):  # this must be modified so that it does a biased above outlier rejection
    p = 50  # PERCENTILE TO USE for split
    perc = np.percentile(data, p)
    upper_half = data[data > perc]
    d = np.abs(data - np.median(upper_half))
    d2 = np.abs(upper_half - np.median(upper_half))
    mdev = np.median(d2)
    s = d / mdev if mdev else 0.
    return s < m


def gauss(x, *p):
    """
        Simple Gaussian function

    :param x: ind var
    :param p: coefficients A, mu, sigma
    :return: numpy array gauss(x)
    """
    A, mu, sigma = p
    return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


def csq_red(model, data, dof=3.):
    """
        Computed a reduced cui square fit of data to model
        Assumes model sums to one

    :param model: expectation value of model
    :param data: observed data
    :param dof: number of degrees of freedom
    :return: Reduced chi square (float)
    """
    total = np.sum(data)
    csq = np.power(np.array(model) * total - np.array(data), 2)
    error = np.array(data)
    error[np.where(error < 9)] = 9
    csq = csq / error
    csq = np.sum(csq)
    csq /= len(data) - dof
    return csq


def minimal_csq(coeffs, data, dof=3., n_it=20, min_thresh=0.005):
    """
        Does a binary search to find a minimal chi square.

    :param coeffs: gaussian fit coefficients
    :param data: column to compute fit
    :param dof: number of DOF for reduced chi square computation
    :param n_it: number of iterations of Binary Search
    :param min_thresh: difference in BS iterations deemed "sufficiently close"
    :return: minimal value of chi square given n_it and min_thresh
    """
    # TODO - anonymize from specific function type
    indices = np.arange(len(data))
    ub = coeffs[1] + 0.9
    lb = coeffs[1] - 1
    ctr = n_it
    csq_past = 100
    csq_now = 0
    quick_csq = lambda x: csq_red(gauss(indices, *[coeffs[0], x, coeffs[2]]), data, dof=dof)
    while ctr and (csq_past - csq_now > min_thresh or ctr == n_it - 1):
        csq_past = csq_now
        midpoint = (ub + lb) / 2.
        l_midpoint = (lb + midpoint) / 2.
        r_midpoint = (ub + midpoint) / 2.
        csq_l = quick_csq(l_midpoint)
        csq_r = quick_csq(r_midpoint)
        if csq_r < csq_l:
            lb = midpoint
            csq_now = csq_r
        else:
            ub = midpoint
            csq_now = csq_l
        ctr -= 1
    midpoint = (ub + lb) / 2.
    return csq_red(gauss(indices, *[coeffs[0], midpoint, coeffs[2]]), data, dof=dof)


def makegaussian(size, fwhm=3, center=None):
    """
    Adapted from Andrew Giessel on GitHub

    Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:, np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4 * np.log(2) * ((x - x0) ** 2 + (y - y0) ** 2) / fwhm ** 2)



def find_nearest(array, value):
    """Finds nearest value in array"""
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def pythag_dist_kernel(size=3):
    s = size + size + 1
    offsetx = np.repeat(np.reshape(range(s), (s, 1)), s, axis=1)
    offsety = np.transpose(offsetx)
    return np.square(offsetx - size) + np.square(offsety - size)


def get_header_info(rolist, info=['RA'], loc='/mir3/iodfitsdb'):
    """
    Get header info for HIRES
    :param rolist: 
    :param info: 
    :param loc: 
    :return: 
    """
    targ = 'Target'
    info.insert(0, targ)
    out = []
    print (info)
    for i in rolist:
        t = fits.open(loc + '/rj' + str(i[0]) + '.' + str(i[1]) + '.fits')[0]
        prt = [t.header['TARGNAME'].strip(' ')]
        for j in info[1:]:
            prt.append(t.header[j])
        print (prt)
        out.append(prt)
    return out


def conv(arr, kernel=[[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]]):
    return convolve2d(arr, kernel)

# ------------------------------------------------------------------------------
#  Utilities for fitting and interpolation
# ------------------------------------------------------------------------------



def spline_interpolate(arr, nseg=20, percentile=50):
    """
    Returns the cubic spline describing data, split into nseg segments. 
    :param arr: input array
    :param nseg: number of same-size segments
    :param percentile: what percentile of each segment to use as midpt y value?
    :return: scipy.interpolate.CubicSpline
    """
    l = len(arr)
    midpts = [np.median(arr[:int(l // nseg // 2)])]
    x = [0]
    for i in range(nseg):
        x += [int((0.5+i)*(l/nseg))]
        midpts += [np.percentile(arr[i * (l // nseg):(i + 1) * (l // nseg)], percentile)]
    x += [l-1]
    midpts += [np.median(arr[-int(l // nseg / 2):])]
    return CubicSpline(x,midpts)

def poly_interpolator(arr,degree=4, nseg=5,percentile=95):
    """
    Returns the function describing a polynomial fitting the data, split into nseg segments
    :param arr: input array
    :param degree: degree of polynpomial fit
    :param nseg: number of segments of equal length to use.
    :param percentile: 
    :return: 
    """
    l = len(arr)
    midpts = [np.median(arr[:int(l // nseg // 2)])]
    x = [0]
    for i in range(nseg):
        x += [int((0.5 + i) * (l // nseg))]
        midpts += [np.percentile(arr[i * (l // nseg):(i + 1) * (l // nseg)], percentile)]
    x += [l - 1]
    midpts += [np.median(arr[-int(l // nseg // 2):])]
    return np.poly1d(np.polyfit(x, midpts, deg=degree))


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')

# This is the most effective continuum fit
def continuum_fit(arr, percentile_kernel = 101,savitzky_kernel = 2001, savitzky_degree=4, perc=50):
    # This fixes the singularities
    # Value of 500 chosen arbitrarily - should not have too much of an effect
    fixval = np.max([np.abs(np.min(arr) * 2),500.])
    fix = arr + fixval
    pcf = percentile_filter(fix, perc, size=percentile_kernel)
    sav = savitzky_golay(pcf, savitzky_kernel, savitzky_degree)
    return fix/(sav/np.max(sav)) - fixval

def deblaze(arr, method = 'savitzky', percentile_kernel = 101, savitzky_kernel=2001, savitzky_degree=4, perc=50):

    if method == 'savitzky':
        return continuum_fit(arr, percentile_kernel=percentile_kernel, savitzky_kernel=savitzky_kernel,
                             savitzky_degree=savitzky_degree, perc=perc)
    elif method == 'meanshift':
        median_of_array = np.median(arr)
        bandwidth = estimate_bandwidth(arr[:, np.newaxis], quantile=0.1)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(arr[:, np.newaxis])
        # xvals = np.arange(4608)
        test = np.array(arr)
        labels = ms.labels_
        # Replace the missing values (not at maximal cluster) with median of array values in cluster in original array
        test[labels != 0] = np.median(test[labels == 0])

        med_test = sg.medfilt(test, kernel_size=101)

        return arr / med_test * median_of_array
    else:
        raise KeyError('The deblaze method you have passed is not implemented. Please pick from savitzky, bstar, and meanshift')



# ------------------------------------------------------------------------------
#  Utilities for laser search
# ------------------------------------------------------------------------------

def getpercentile(order, perc, method='meanshift', kernel_bandwidth=100, kernel='epanechnikov'):
    """
    Returns value of 'perc'th percentile
    (usually 75th) count value in 'order'

    :param order: Spectral order to compute percentile on
    :param perc: What(th) %ile to compute.
    """
    #TODO - add support for piecewise thresholds

    if method == 'percentile':

        nsplits = 1  # Compute percentile piecewise - I have not been
        maximum_thresh = 0
        l=len(order)
        inc = l / nsplits
        for i in range(nsplits):
            sub = order[i * inc:(i + 1) * inc]
            percentile = np.percentile(sub, perc)
            if maximum_thresh < percentile:
                maximum_thresh = percentile
        return maximum_thresh
    elif method == 'kde':
        kde = KernelDensity(kernel=kernel, bandwidth=kernel_bandwidth).fit(order)

    elif method == 'meanshift':
        bandwidth = estimate_bandwidth(order[:,np.newaxis], quantile=0.1)
        # print ('Bandwidth is {0}'.format(bandwidth))
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(order[:,np.newaxis])
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)
        #for k in range(n_clusters_):
        #    my_members = labels == k
        #   print "cluster {0}: {1}".format(k, order[:, np.newaxis][my_members, 0])
        #return cluster_centers[0][0]
        # print(cluster_centers)
        # TODO Determine why there is a index error ocurring here - should there be more than 3 clusters
        # or is this normal behavior?
        try:
            return np.max([cluster_centers[0][0],cluster_centers[1][0],cluster_centers[2][0]])
        except IndexError:
            try:
                return np.max([cluster_centers[0][0], cluster_centers[1][0]])
            except IndexError:
                return cluster_centers[0][0]
    else:
        raise KeyError




def contiguous_regions(condition):
    """
    Borrowed from Joe Kington on StackOverflow

    Finds contiguous True regions of the boolean array 'condition'. Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.
    """

    # Find the indices of changes in 'condition'
    d = np.diff(condition)
    idx, = d.nonzero()

    # We need to start things after the change in 'condition'. Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]  # Edit

    # Reshape the result into two columns
    idx.shape = (-1, 2)
    return idx


def finddeviates(order, thresh, npix=3):
    """returns a list of deviation indices [start,stop]."""
    out = []
    #plt.plot(order)
    #plt.show()
    #print(order, thresh)
    #plt.cla()
    for start, stop in contiguous_regions(order > thresh):
        diff = stop - start
        if diff >= npix:
            out.append([diff, start])
    return out


def findthresh(order, npix=3.0, method='full_deviation'):
    """
    Computes threshold using median absolute deviation (MAD) method

    note here that 'order' passed to this function is generally:
    spectral order - 75th percentile count value of spectral order

    That is, the sum of 3 consecutive points in order being above 0
     is the same as sum of 3 being above 75th %ile in spectral order.

    Returns the Median Absolute (positive) Deviation of npix (usually 3)
     pixel bins above the percentile set in findhigher.
    """
    # Number of pixels to demand consecutively deviant. 3 is appropriate for HIRES.
    # Convolution as a way of binning by 3 pixels to see groups that exceed
    if method == '3pix':
        binned_ord = np.convolve(np.ones(npix) / npix, order, 'valid')
        deviations = finddeviates(binned_ord, 0, npix)
        uppies = []
        for j in deviations:
            for i in range(j[0]):
                uppies.append(binned_ord[j[1] + i])
    elif method == 'full_deviation':
        deviations = finddeviates(order, 0, npix)
        uppies = []
        for j in deviations:
            uppies.append(np.median(order[j[1]:j[1] + j[0]]))
    return np.median(uppies)




def has_singularity(order):
    """ Tests order for a 'Singularity', or a common reduction error resulting in huge counts"""
    order = order[4:-4]
    order_c = np.convolve(np.abs(order), [1, 1])
    big = np.where(order_c > 500)[0]
    zero_crossings = np.where(np.diff(np.sign(order)))[0]
    return True if np.intersect1d(big, zero_crossings).__len__() else False


def hires_ignored_wavelengths(wav):
    """ Ignore results from these wavelength ranges."""
    ignore_lb = np.array([6557.78, 4855.78, 4335.78, 3964.78, 3928.46, 3963.25,
                          5890.7, 5884.73, 7585.0, 7964.0, 6863.0])
    ignore_ub = np.array([6568.22, 4866.22, 4346.22, 3975.22, 3938.9, 3973.69,
                          5901.14, 5895.17, 7660.0, 7965.2, 6920.0])
    onesarr = np.ones(len(ignore_lb))
    return np.any(np.logical_xor(np.greater(ignore_lb, wav * onesarr), np.greater(ignore_ub, wav * onesarr)))


#  TODO ORGANIZE, should be in Output


def view_raw(raw):
    # Takes a raw spectrum and plots the entire CCD image
    plt.imshow(raw.data,vmax=np.percentile(raw.data,99),vmin=np.percentile(raw.data,2),interpolation='nearest')



# from Eelco Hoogendoorn on stackoverflow
# https://stackoverflow.com/questions/21030391/how-to-normalize-an-array-in-numpy
def normalized(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)

def median_of_one(arr):
    med = np.median(arr)
    return arr/med
