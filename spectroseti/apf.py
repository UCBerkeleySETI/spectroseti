# apf.py
#
# Nate Tellis 2017
#
#
# Extends class definitions in spectra.py for the High Resolution Echelle
#    spectrometer at the WM Keck observatory.

__author__ = 'nate'
import numpy as np
from astropy.io import fits
from os import listdir
from astropy.stats import median_absolute_deviation as mad_std
from scipy.stats import mode
import random
import scipy.ndimage.filters as filters
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

import apfdefinitions as apfdefs
import definitions as defs
import utilities as utilities
import spectra as spec
from tqdm import tqdm



class APFRedObs(spec.ReducedObs):

    raw_red_correspondence = None

    def __init__(self, run, obs, atlas=spec.WavelengthAtlas()):
        self.devs_set = 0
        self.run = ''
        self.obs = 0
        self.loadobservation(run, obs)
        self.atlas = atlas

    def loadobservation(self, run, obs, deblaze=0):
        """
        Load a reduced APF spectrum into object params
        :param run: 3 character string, i.e. 'aii', denoting observing run
        :param obs: int, denoting observation within a run
        :return: none
        """

        self.devs_set = 0
        self.run = run
        self.obs = obs

        try:
            self.dat = fits.open(apfdefs.spectra_dir_reduced+'r%(run)s.%(obs)s.fits' % locals())
            self.wav = fits.open(defs.project_root_dir+apfdefs.apf_wavs)
            self.wavs = self.wav[0].data
            self.counts = self.dat[0].data
            if deblaze:
                self.deblaze_orders('savitzky')

        except IOError:
            self.dat = None
            print('Fatal Error - r%(run)s.%(obs)s.fits not found!' % locals())

    def deblaze_orders(self,method='meanshift', percentile_kernel = 101, savitzky_kernel=2001,
                       savitzky_degree=4, perc=50, bstar_correction=None):
        if method == 'savitzky':
            deblaze = lambda x: utilities.deblaze(x, method='savitzky', percentile_kernel=percentile_kernel,
                                                  savitzky_kernel=savitzky_kernel,
                                                  savitzky_degree=savitzky_degree, perc=perc)
            self.counts = np.apply_along_axis(deblaze, 1, self.counts)
        elif method == 'bstar':
            try:
                bstar_correction.shape
            except AttributeError:
                bstar_correction = np.load(defs.project_root_dir+apfdefs.bstar_correction_dir)
            self.counts = self.counts / bstar_correction

        # Meanshift deblaze only really works once another deblaze has been applied
        elif method == 'meanshift':
            deblaze = lambda x: utilities.deblaze(x, method='meanshift', percentile_kernel=percentile_kernel,
                                                  savitzky_kernel=savitzky_kernel,
                                                  savitzky_degree=savitzky_degree, perc=perc)
            self.counts = np.apply_along_axis(deblaze, 1, self.counts)
        else:
            raise KeyError(
                'The deblaze method you have passed is not implemented. Please pick from savitzky, bstar, and meanshift')



    def loaddevs(self,method='simple',n_mads=5,percentile=75):
        if method=='simple':
            self.devs, self.percentiles_and_thresholds = findhigher(self, n_mads, percentile, atlas=self.atlas)
            self.devs_set = 1
        # An even simpler method that is only 4 pix >1.3*med
        if method=='simpler':
            self.devs, self.percentiles_and_thresholds  = findhigher(self, n_mads, percentile, atlas=self.atlas, method='basic')
            self.devs_set = 1
        if method == 'MAD':
            raise NotImplementedError

    def red2raw(self, ord, pix):

        try:
            if not self.raw_red_correspondence:
                self.raw_red_correspondence = apfdefs.correspondence
        except IOError:
            print('Fatal Error - apf order correspondence mask not found!')

        return self.raw_red_correspondence[ord,pix]


class APFRawObs(spec.RawObs):
    """
    Extends RawObs class for the specifics of the APF
    """
    data = None
    raw_red_correspondence = None
    rawdims=(4608,2080)

    def __init__(self, run, obs):
        self.load(run,obs)
        self.raw_red_correspondence = np.load(defs.project_root_dir + apfdefs.correspondence)

    def load(self, run, obs):
        """Load a raw observation."""
        raw_fits = fits.open(apfdefs.spectra_dir_raw + 'ucb-%(run)s%(obs)s.fits' % locals())
        self.data = np.rot90(raw_fits[0].data,3) # rotation to get raw in line with prefered viewing orientation
        raw_fits.close()

    def show_raw(self):
        fig, ax = plt.subplots()
        ax.imshow(np.transpose(raw[0].data), vmin=np.min(raw[0].data), vmax=np.percentile(raw[0].data, 85))

    def red2raw(self, ord, pix):
        '''try:
            if self.raw_red_correspondence == None:
                self.raw_red_correspondence = np.load(defs.project_root_dir + apfdefs.correspondence)
        except IOError:
            print('Fatal Error - apf order correspondence mask not found!')'''

        return self.raw_red_correspondence[ord, pix]

    def retrieve_subset(self, ord, index, yradius=5, xradius=5):
        """
        Retrieves a rectangular subset of the raw image at
        the stellar ridge corresponding to some order/pixel pair
        """
        central = self.red2raw(ord,index)
        ylow = 0 if central - yradius < 0 else central - yradius
        yhigh = self.rawdims[1] if central + yradius > self.rawdims[1] else central + yradius
        xlow = 0 if index - xradius < 0 else index - xradius
        xhigh = self.rawdims[0] if index - xradius > self.rawdims[0] else index + xradius
        return self.data[int(ylow):int(yhigh),int(xlow):int(xhigh)]



    def cr_reject(self, order, pix):
        raise NotImplementedError

    def run_obs_to_filename(self, run, obs):
        """Simple utility to translate from run/obs int pair to raw filename."""
        return 'ucb-%(run)s%(obs)s.fits' % locals()



# DEPRECATED
'''def populateapfobsdatabase(filename='apf_obs.db'):
    filenames = listdir(apfdefs.logsheet_dir)
    log_accumulator = []
    for log in filenames: # Replace with all logsheets
        run = log.split('.')[0]
        f = open(apfdefs.logsheet_dir+log)
        logs = f.readlines()
        for i in range(7,13):
            try:
                if logs[i][:3].isdigit():
                    break
            except IndexError:
                continue
        log_accumulator += list(map(lambda x: str.split(run+' '+x,maxsplit=12), logs[i:]))
        f.close()
        length = len(sorted(log_accumulator,key=len, reverse=True)[0])
        acc_array = np.array([xi+[None]*(length-len(xi)) for xi in log_accumulator])'''


# Turned off AltMinThresh for now. See what this does.
def find_deviations(ords, wavs, order, perc=75, n_mads=5,alt_min_thresh=0, atlas=spec.WavelengthAtlas(), out=[], npix=3, acc=[]):
    """
        The meat of the laser search pipeline, this function finds all pixel regions that deviate
        over an entire spectroscopic order.

    :param ords: Order counts , numpy.ndarray
    :param wavs: Corresponding wavelength range, numpy.ndarray
    :param order: Which order to find deviations in?
    :param simulator: DEPRECATED - for adding simulated lasers to a run
    :param perc: What(th) percentile should the MAD be computed above? For continuum determination
    :param n_mads: How many MADs should the code look for deviations above?
    :param atlas: Wavelength atlas for night sky and emission line removal
    :param out: Can pass a list for accumulation
    :param npix: How many pixels to demand consecutive
    :return: list of lists with format [[order,start,len,nth_percentile, threshold set, mean deviate pixel
                                        median deviate pixel, midpt pixel value], ...]
    """
    l = len(ords[0])
    percentile = utilities.getpercentile(ords[order], perc)
    threshold = utilities.findthresh(ords[order] - percentile)
    th2 = 100 + percentile  # TODO - fix second thresholding
    if percentile < 100 and threshold < th2:
        final_threshold = percentile + n_mads * threshold * 1.5
        secondary_threshold = percentile + n_mads * threshold / 2.0
    else:
        final_threshold = percentile + n_mads * threshold
        secondary_threshold = percentile + n_mads * threshold / 2.0
    # Experimental basic thresholding so that MADS is not far too high
    if percentile > 500 and n_mads * threshold > 0.3 * percentile and alt_min_thresh:
        final_threshold = 1.3 * percentile
        secondary_threshold = 1.15 * percentile
    # Accumulator for testing whole-dataset thresholding
    acc.append([percentile, final_threshold])
    contig = utilities.finddeviates(ords[order], final_threshold, npix=npix)
    if len(contig) != 0:
        for x in contig:
            midpt = x[0] // 2 + x[1]
            if l*0.02 < midpt < l-l*0.02 and x[1]>l*0.01 and x[0]+x[1]<l-l*.01: #\ COMMENTED OUT THE IGNORED WAVELENGTHS
                   # and not hires_ignored_wavelengths(wavs[order][midpt]) \
                   # and not list(atlas.ns_lookup(wavs[order][midpt])) and not has_singularity(ords[order]):
                deviation_pixels = ords[order][x[1]:x[1] + x[0]] - final_threshold
                out.append(
                    [order, x[0], x[1], float(percentile), float(threshold), float(np.mean(deviation_pixels)),
                     float(np.median(deviation_pixels)), float(wavs[order][midpt])])
    return out

def find_deviations_basic(ords, wavs, order, perc=75, n_mads=5,alt_min_thresh=1, atlas=spec.WavelengthAtlas(), out=[], npix=3, acc=[]):
    """
        The meat of the laser search pipeline, this function finds all pixel regions that deviate
        over an entire spectroscopic order.

    :param ords: Order counts , numpy.ndarray
    :param wavs: Corresponding wavelength range, numpy.ndarray
    :param order: Which order to find deviations in?
    :param simulator: DEPRECATED - for adding simulated lasers to a run
    :param perc: What(th) percentile should the MAD be computed above? For continuum determination
    :param n_mads: How many MADs should the code look for deviations above?
    :param atlas: Wavelength atlas for night sky and emission line removal
    :param out: Can pass a list for accumulation
    :param npix: How many pixels to demand consecutive
    :return: list of lists with format [[order,start,len,nth_percentile, threshold set, mean deviate pixel
                                        median deviate pixel, midpt pixel value], ...]
    """
    l = len(ords[0])
    percentile = utilities.getpercentile(ords[order], perc)
    threshold = utilities.findthresh(ords[order] - percentile)
    final_threshold = 1.5 * percentile
    # Accumulator for testing whole-dataset thresholding
    if order==20:
        pass
    acc.append([percentile, final_threshold])
    contig = utilities.finddeviates(ords[order], final_threshold, npix=npix)
    if len(contig) != 0:
        for x in contig:
            midpt = x[0] // 2 + x[1]
            if x[1] > 500 and x[0]+x[1] < 4108 and x[0] < 10: #\ COMMENTED OUT THE IGNORED WAVELENGTHS
                   # and not hires_ignored_wavelengths(wavs[order][midpt]) \
                   # and not list(atlas.ns_lookup(wavs[order][midpt])) and not has_singularity(ords[order]):
                deviation_pixels = ords[order][x[1]:x[1] + x[0]] - final_threshold
                out.append(
                    [order, x[0], x[1], float(percentile), float(threshold), float(np.mean(deviation_pixels)),
                     float(np.median(deviation_pixels)), float(wavs[order][midpt])])
    return out


def findhigher(obs, n_mads, perc, atlas=spec.WavelengthAtlas(),method='original'):
    """
    Finds potential laser candidates using Median
    Absolute Deviation method in find_deviations,
    potentially adding simulated signals as well.

    First bins spectrum, having subtracted 75th percentile
    then finds median 4-pix deviation, uses this as threshold
    by which to search for every multi-pixel deviation
    of at least 4 pixels


    :param atlas: Wavelength atlas for throwing out night sky line frequencies
    :param obs: HiresObs instance
    :param n_mads: Number of MADs above %ile to use as threshold
    :param perc: Percentile to use for stellar continuum
    :return: list of lists with format [[order,start,len,nth_percentile, threshold set, mean deviate pixel
                                        median deviate pixel, midpt pixel value], ...]
    """
    out = []
    per_thr_accumulator = []
    cts = obs.counts
    wavs = obs.wavs
    if method=='original':
        print('Searching order-by-order...')
        for i in tqdm(range(79), leave=False, miniters=8):
            out = find_deviations(cts, wavs, i, perc=perc, n_mads=n_mads,
                                  atlas=atlas, out=out, acc=per_thr_accumulator)
        return out, per_thr_accumulator
    elif method=='basic':
        print('Searching order-by-order...')
        for i in tqdm(range(79), leave=False, miniters=8):
            out = find_deviations_basic(cts, wavs, i, perc=perc, n_mads=n_mads,
                                  atlas=atlas, out=out, acc=per_thr_accumulator,npix=4)
        return out, per_thr_accumulator
    else:
        raise NameError('Incorrect method keyword')


def deblaze(arr):
    l = len(arr[:-1])
    CS = utilities.poly_interpolator(arr[:-1],degree=5,nseg=15,percentile=80)
    return arr[:-1]/(CS(np.arange(l))/np.percentile(CS(np.arange(l)),95))


