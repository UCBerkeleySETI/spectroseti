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




class SearchOutputManager:


    def __init__(self):
