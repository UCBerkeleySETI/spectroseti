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
import sqlite3


# Class manages IO for laser search
#
# writes/reads output from searches, reads logsheets, maintains stats DB,

class LaserSearchIOManager:

    def __init__(self):
        pass


    # Write searchlog



    # Read searchlog



    # Find (run,obs) pairs matching:
    # - date or daterange
    # - targname or array of targnames
    # - run or runrange


    # Check existence of (run.obs) pair
    # TODO - this should be a method in utils, probably, generalized for {fileexists}