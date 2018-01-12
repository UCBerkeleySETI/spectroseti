# apf.py
#
# Nate Tellis 2017
#
#
# Extends class definitions in spectra.py for the High Resolution Echelle
#    spectrometer at the WM Keck observatory.

__author__ = 'nate'


import definitions as defs
import utilities as utilities
import spectra as spec
import apf as apf
import apfdefinitions as apfdefs
import output


# highest level laser search class. Should contain helper methods to extract specific targets etc
class LaserSearch():

    raw_directory = apfdefs.spectra_dir_raw
    reduced_directory = apfdefs.spectra_dir_reduced

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
        reduced_spectrum.loaddevs(method=load_devs_method,n_mads=number_mads,percentile=search_percentile)
        return reduced_spectrum


    def search_multiple(self, observations, output_pngs=0, logfile=0, db_write=0, stats=0):
        # observations expects a tuple of run,obs pairs

        # setup directories, filenames, local accumulator variables etc
        if output_pngs or logfile:
            pass

        # a little pseudocodey
        for observation in observations:
            run = observation[0]
            obs = observation[1]
            reduced_spectrum = self.search_one(run, obs)

            if output_pngs:
                # generate and save all output

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