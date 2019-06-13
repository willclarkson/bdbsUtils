#
# txt2fits.py 
#

# WIC 2019-06-12 - steps to ingest bdbs v2 into fits files, with
# HEALPIX-denoted object IDs

import os, time
import numpy as np
import healpy

class bdbsCat(object):

    """Object for BDBS catalog"""

    def __init__(self, infil='TEST.catalog', nsidePow=22):

        self.infil=infil[:]

        #catalog array
        self.aCat=np.array([])

        # power of 2 for NSIDE parameter
        self.nsidePow = np.copy(nsidePow)
        self.setNSIDE()

        # RA, DEC columns
        self.colRA = 0
        self.colDE = 1

        self.ra = np.array([])
        self.de = np.array([])

        # ID column
        self.id = np.array([])

        # Control variable
        self.Verbose=True

    def loadCat(self):

        """Loads the catalog"""

        if not os.access(self.infil, os.R_OK):
            print("bdbsCat.loadCat FATAL - cannot read path: %s" \
                  % (self.infil))
            return

        self.aCat = np.loadtxt(self.infil, unpack=False)
        
    def setCOO(self):

        """Sets the sky coordinates from the input data"""
        
        if np.size(self.aCat) < 1:
            if self.Verbose:
                print("bdbsCat.setCOO FATAL - catalog not yet imported")
            return

        self.ra = self.aCat[:,self.colRA]
        self.de = self.aCat[:,self.colDE]

    def estHealpixSpacing(self):

        """Rough estimate of healpix spacing from nside power"""

        estSpacing = 58.6 * 3600. * 2.0**(0.-self.nsidePow)

        print("bdbsCat.estHealpixSpacing INFO - at nsidepow %i, NSIDE=%i, mean spacing roughly %.2e arcsec (%.2e degrees)" % (self.nsidePow, 2.0**(self.nsidePow), estSpacing, estSpacing/3600.))

    def setNSIDE(self):

        """Sets NSIDE from the nsidepower"""

        self.nside = long(2.0**(self.nsidePow))

    def setID(self):

        """Creates ID from the coordinates"""

        self.id = healpy.ang2pix(self.nside, self.ra, self.de, lonlat=True)
        
    def wrapLoadAndID(self):

        """Wrapper - loads the catalog and sets the ID"""

        self.setNSIDE()
        self.loadCat()
        self.setCOO()
        self.setID()

def go(nsidePow=10):

    """Tester for some healpix functionality. Put nsidePow<=10 to reproduce the NASA HEALPIX frontpage documenation examples."""

    BD = bdbsCat(nsidePow=nsidePow)
    BD.estHealpixSpacing()
    BD.wrapLoadAndID()

