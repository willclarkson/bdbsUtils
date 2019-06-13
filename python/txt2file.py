#
# txt2fits.py 
#

# WIC 2019-06-12 - steps to ingest bdbs v2 into intermediate files
# with column names for later ingestion into database system. Object
# IDs are assigned via HEALPIX.


import os, time
import numpy as np
import healpy

class bdbsCat(object):

    """Object for BDBS catalog"""

    def __init__(self, infil='TEST.catalog', nsidePow=22):

        self.infil=infil[:]

        # output catalog name
        self.outfil='TEST.csv'
        self.outExt='csv'
        self.outDir='tiles'

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

    def loadEntireCat(self):

        """Loads the catalog"""

        if not os.access(self.infil, os.R_OK):
            print("bdbsCat.loadEntireCat FATAL - cannot read path: %s" \
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

        """Wrapper - loads the entire catalog and sets the ID"""

        self.setNSIDE()
        self.loadEntireCat()
        self.setCOO()
        self.setID()

    def setOutputFilename(self, nMin=3):

        """Sets the output filename based on the input filename"""
        
        if len(self.infil) < nMin:
            if self.Verbose:
                print("bdbsCat.setOutputFilename FATAL - refusing to operate on filename shorter than %i characters" % (nMin))
            return

        self.outfil="%s.%s" % (os.path.splitext(self.infil)[0], self.outExt)


    def ensureOutdirExists(self):

        """Ensures the output directory exists"""

        if not os.path.isdir(self.outDir):
            os.makedirs(self.outDir)

    def processStream(self, iMax=-1, DBG=False):

        """Reads the input file line by line, operating on what it finds"""

        if not os.access(self.infil, os.R_OK):
            print("bdbsCat.processStream FATAL - cannot read path %s" % (self.infil))
            return

        if len(self.outfil) < 2:
            print("bdbsCat.processStream FATAL - output file not set.")
            return
        
        pathOut="%s/%s" % (self.outDir, self.outfil)

        iCount = 0
        with open(self.infil, "r") as rObj, open(pathOut, "w") as wObj:
            for line in rObj:

                lProc = self.processInputLine(line, DBG=DBG)

                wObj.write("%s" % (lProc))

                # counter so that we only need open the first few
                # lines for testing
                iCount = iCount + 1
                if iCount > iMax and iMax > 0:
                    break

    def processInputLine(self, lineIn='', DBG=False):

        """For a single input string, extracts the coordinates,
        assigns a HEALPIX-based ID, tacks the ID onto the beginning of
        the output string, and returns the string including
        commas. Note that the carriage return is retained on the
        output as well as the input."""

        lineOut=lineIn[:]

        if len(lineOut) < 1:
            return lineIn

        lSplit = lineIn.split()
        
        # nothing to be done if it's impossible to have coordinates
        # here.
        if len(lSplit) < 2:  
            return lineIn

        # extract the coordinates and assign a HEALPIX ID
        thisRA = np.float(lSplit[self.colRA])
        thisDE = np.float(lSplit[self.colDE])

        # hack since I don't have updated healpy on my laptop
        try:
            thisID = healpy.ang2pix(self.nside, thisRA, thisDE, lonlat=True)
        except:
            thisID = healpy.ang2pix(self.nside, np.radians(thisRA)-np.pi, np.radians(thisDE))

        lSplit.insert(0,"%i" % (thisID))

        lineOut=','.join(lSplit)
        lineOut='%s\n' % (lineOut)

        if DBG:
            print lSplit
            print lineOut

        return lineOut

def go(nsidePow=10):

    """Tester for some healpix functionality. Put nsidePow<=10 to reproduce the NASA HEALPIX frontpage documenation examples."""

    BD = bdbsCat(nsidePow=nsidePow)
    BD.estHealpixSpacing()
    BD.wrapLoadAndID()

def testStream(nMax=5):

    """Tests stream-based reading of input file"""
    
    BD = bdbsCat(nsidePow=22)
    BD.setOutputFilename()
    BD.ensureOutdirExists()

    print BD.outfil
    BD.processStream(iMax=nMax, DBG=True)

    # this works reasonably well. To add:
    #
    # (i) bunching of output lines to avoid HDD thrashing
    #
    # (ii) column names entry for csv files
    # 
    # (iii) Option to set max lines for output files
