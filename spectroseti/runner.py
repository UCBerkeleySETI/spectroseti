# apf.py
#
# Nate Tellis 2017
#
#
# Extends class definitions in spectra.py for the High Resolution Echelle
#    spectrometer at the WM Keck observatory.
import spectroseti.output

__author__ = 'nate'


import definitions as defs
import utilities as util
import spectra as spec
import apf as apf
import apfdefinitions as apfdefs
#import output
import copy
import numpy as np
from tqdm import tqdm


# highest level laser search class. Should contain helper methods to extract specific targets etc
class LaserSearch():

    raw_directory = apfdefs.spectra_dir_raw
    reduced_directory = apfdefs.spectra_dir_reduced
    bstar_correction = np.load(defs.project_root_dir + apfdefs.bstar_correction_dir)

    # Default initialization
    # sets raw, red directories to those in apf_def
    def __init__(self, raw_dir = apfdefs.spectra_dir_raw, red_dir = apfdefs.spectra_dir_reduced):
        self.raw_directory = raw_dir
        self.reduced_directory = red_dir
        pass


    def search_one(self, run, observation, load_devs_method='simple',number_mads=5,search_percentile=75):
        # don't need the raw at first
        # raw = apf.APFRawObs(run, observation)
        reduced_spectrum = apf.APFRedObs(run, observation)

        # Now first deblaze (Savizky works better on lower orders than b-star)
        #reduced_spectrum.deblaze_orders(method='bstar',bstar_correction=self.bstar_correction)
        reduced_spectrum.deblaze_orders(method='savitzky')

        # Make a copy of the (bstar-only) deblazed spectrum
        #bstar_deblazed = copy.deepcopy(reduced_spectrum)
        print('Beginning Meanshift deblazing')
        # Meanshift deblaze the reduced spectrum
        reduced_spectrum.deblaze_orders(method='meanshift')
        # Load deviations with the meanshift method
        # loaddevs -> findhigher -> find_deviations -> getpercentile (has meanshift method)
        reduced_spectrum.loaddevs(method=load_devs_method,n_mads=number_mads,percentile=search_percentile)

        # Here we go back and check the bstar spectrum for the same positives
        # One way to proceed is:
        #   compute perc from bstar_deblazed
        #   compute thresh
        #   ensure three pixels are > perc_b+n_mads*thresh_b

        return reduced_spectrum


    def search_multiple(self, observations, output_pngs=0, logfile=0, db_write=0, stats=0):
        # observations expects a tuple of run,obs pairs
        # setup directories, filenames, local accumulator variables etc
        if stats:
            positive_accumulator = np.zeros(79)

        if output_pngs or logfile:


            pass

        # a little pseudocodey
        for observation in observations:
            run = observation[0]
            obs = observation[1]
            reduced_spectrum = self.search_one(run, obs)

            if output_pngs:
                # generate and save all output
                raw = None
                try:
                    raw = apf.APFRawObs(run,obs)
                except:
                    raw = None
                ndev = len(reduced_spectrum.devs)
                print('Writing output images to '+ apfdefs.output_png_dir)
                for i in tqdm(range(ndev), miniters=int(ndev/10)):
                    # TODO pass down a folder here for saving the output
                    spectroseti.output.view_dev(reduced_spectrum, devnum=i, raw=raw, save=1)
                pass

            #TODO this is first priority
            if logfile:
                # Accumulate statistics for logfile
                pass

            if db_write:
                # Save to database
                pass

        # write logfile
        if logfile:
            #save logfile to logfile directory, as well as search run directory
            pass