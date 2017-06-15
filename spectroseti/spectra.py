# spectra.py
#
# Nate Tellis 2017
#
#
# Contains template classes and methods for manipulating spectra

import numpy as np
import definitions as defs
import utilities as utilities
from utilities import getpercentile, findthresh, finddeviates, hires_ignored_wavelengths, has_singularity

__author__ = 'nate'

from astropy.io import fits
from astropy.stats import median_absolute_deviation as mad_std
from scipy.stats import mode
import random
import scipy.ndimage.filters as filters
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt




# ------------------------------------------------------------------------------
#  Class definitions
# ------------------------------------------------------------------------------


class WavelengthAtlas:
    """
    Simple class for loading and holding the night sky line wavelength atlas
    as well as any other wavelengths that need to be maintained.
    """

    def __init__(self):
        self.load_night_sky()

    def load_night_sky(self, filename='default'):
        """Loads night sky atlas datafile"""
        if filename == 'default':
            self.ns_atlas = np.load(defs.project_root_dir + defs.ns_atlas)
        else:
            self.ns_atlas = np.load(filename)

    def ns_lookup(self, wavelength, tolerance=0.30):
        """Lookup set of test wavelengths closes with a given tolerance, see if night sky"""
        # TODO - Possible to make this check if brightness of night sky is expected
        near = utilities.find_nearest(self.ns_atlas[:, 0], wavelength)
        if abs(near - wavelength) < tolerance:
            return self.ns_atlas[self.ns_atlas[:, 0] == near][0]
        else:
            return []


class ReducedObs:
    """
    Template container class for reduced echelle spectra.
    Assumes a run/obs indexing framework
    
    """
    
    def __init__(self, run, obs, atlas=WavelengthAtlas()):
        self.devs_set = 0   # Has a laser search been run on this observation?
        self.run = 0
        self.obs = 0
        self.loadobservation(run, obs)
        self.atlas = atlas

    def loadobservation(self, run, obs):
        raise NotImplementedError

    def loaddevs(self):
        raise NotImplementedError


class RawObs:
    """
    Class encapsulates raw spectrum
    and associated metadata. Contains relevant
    routines for manipulation of raw data.

    Uses a file 'yinds' containing a mapping from
    (order, pixel) --> (row, column)
    on the raw image, such that the resulting
    indices lie on the stellar ridge.
    
    For use with a new spectral survey, extend all methods
    """

    def __init__(self, run, obs, yinds='null'):
        raise NotImplementedError

    def load(self, run, obs):
        """Load a raw observation."""
        raise NotImplementedError

    def retrieve_subset(self, ord, index, yextent=5, xextent=5):
        """
        Retrieves a rectangular subset of the raw image at
        the stellar ridge corresponding to some order/pixel pair
        """
        raise NotImplementedError

    def cr_reject(self, order, pix):
        raise NotImplementedError

    def run_obs_to_filename(self, run, obs):
        """Simple utility to translate from run/obs int pair to raw filename."""
        ostr = str(obs)
        return 'j' + str(run) + np.abs(len(ostr) - 4) * '0' + str(obs) + '.fits'



class LaserHunter:

    def __init__(self,alg=''):
        buffer = (20,20)  # Number of pixels at each end of an order to ignore
        alg=''
        raise NotImplementedError


    def find_deviations(ords, wavs, order, perc=75, n_mads=5,alt_min_thresh=1, atlas=WavelengthAtlas(), out=[], npix=3, acc=[]):
        #   NOTE - CURRENTLY WORKS FOR HIRES

        """
            The meat of the laser search pipeline, this function finds all pixel regions that deviate
            over an entire spectroscopic order.
    
        :param ords: Order counts , numpy.ndarray
        :param wavs: Corresponding wavelength range, numpy.ndarray
        :param order: Which order to find deviations in?
        :param perc: What(th) percentile should the MAD be computed above? For continuum determination
        :param n_mads: How many MADs should the code look for deviations above?
        :param atlas: Wavelength atlas for night sky and emission line removal
        :param out: Can pass a list for accumulation
        :param npix: How many pixels to demand consecutive
        :return: list of lists with format [[order,start,len,nth_percentile, threshold set, mean deviate pixel
                                            median deviate pixel, midpt pixel value], ...]
        """
        l = len(ords[0])
        percentile = getpercentile(ords[order], perc)
        threshold = findthresh(ords[order] - percentile)
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
        contig = finddeviates(ords[order], final_threshold, npix=npix)
        if len(contig) != 0:
            for x in contig:
                midpt = x[0] / 2 + x[1]
                if l*0.02 < midpt < l-l*0.02 and x[1]>l*0.01 and x[0]+x[1]<l-l*.01 and not hires_ignored_wavelengths(wavs[order][midpt]) \
                        and not list(atlas.ns_lookup(wavs[order][midpt])) and not has_singularity(ords[order]):
                    deviation_pixels = ords[order][x[1]:x[1] + x[0]] - final_threshold
                    out.append(
                        [order, x[0], x[1], float(percentile), float(threshold), float(np.mean(deviation_pixels)),
                         float(np.median(deviation_pixels)), float(wavs[order][midpt])])
        return out