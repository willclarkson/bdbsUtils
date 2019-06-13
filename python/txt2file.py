#
# txt2fits.py 
#

# WIC 2019-06-12 - steps to ingest bdbs v2 into intermediate files
# with column names for later ingestion into database system. Object
# IDs are assigned via HEALPIX.


import os, time, sys
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

        self.lineHeader=''

        # are we splitting the output across multiple files?
        self.splitFiles=False

        # number of lines to write in a bunch
        self.nBunch=1000000

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

        # this is no longer required, since the output file is always
        # set below.
        #if len(self.outfil) < 2:
        #    print("bdbsCat.processStream FATAL - output file not set.")
        #    return
        
        # We want this to operate sensibly whether the outputs are to
        # be split across multiple files or not. The easiest way is
        # probably to use the same variable convention for both
        # cases. So:
        iBunch = 0
        pathOut=self.getIthOutfile(iBunch)

        # we write the header line to the file, then later will open it as APPEND
        if len(self.lineHeader) > 1:
            with open(pathOut, "w") as wObj:
                wObj.write("%s\n" % (self.lineHeader))

        # reference time for reporting the timing later on
        t0 = time.time()

        iCount = 0
        with open(self.infil, "r") as rObj:

            # do this in pieces to avoid excessive writeout thrashing
            bunch=[]

            for line in rObj:
                lProc = self.processInputLine(line, DBG=DBG)
                bunch.append(lProc)
                iCount = iCount + 1

                if len(bunch) == self.nBunch:
                    with open(pathOut, "a") as wObj:
                        wObj.writelines(bunch)
                    bunch = []
                    iBunch = iBunch + 1

                    # report progress to screen?
                    if self.Verbose:
                        sys.stdout.write("\r bdbsCat.processStream INFO - written line %-11i of %s to %s after %.2e seconds" \
                                       % (iCount, self.infil, pathOut, time.time()-t0))
                        sys.stdout.flush()
                    # if we're splitting the output up by bunch,
                    # create the new pathname and initialise the file
                    if self.splitFiles:
                        pathOut = self.getIthOutfile(iBunch)
                        with open(pathOut, "w") as wObj:
                            wObj.write("%s\n" % (self.lineHeader))

                # wObj.write("%s" % (lProc))

                # counter so that we only need open the first few
                # lines for testing
                if iCount > iMax and iMax > 0:
                    break

            # once broken out of the loop, ensure the remaining
            # partial bunch is also written to disk.
            if len(bunch) > 0:
                with open(pathOut, "a") as wObj:
                    wObj.writelines(bunch)

            if self.Verbose:
                sys.stdout.write("\r bdbsCat.processStream INFO - written line %-11i of %s to %s after %.2e seconds" \
                                     % (iCount, self.infil, pathOut, time.time()-t0))
                sys.stdout.flush()

    def getIthOutfile(self, iFiles=0, exten='csv'):

        """One-liner to return the current output filename"""

        stem = os.path.splitext(self.infil)[0]
        count = str(iFiles).zfill(3)
        return "%s/%s_%s.%s" % (self.outDir, stem, count, exten)

    def processInputLine(self, lineIn='', DBG=False):

        """For a single input string, extracts the coordinates,
        assigns a HEALPIX-based ID, tacks the ID onto the beginning of
        the output string, and returns the string including
        commas. Note that the carriage return is retained on the
        output as well as the input.

        At the moment this is pretty dumb, and may return unusual line
        lengths e.g. if the ID setting process fails. For now, we
        trust the input data to be uniform."""

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

    def setHeaderLineCSV(self):

        """Sets header line for CSV. Currently this is hardcoded
        because I only have the one application in mind."""

        listHeader=['id', 'ra', 'dec', 'raerr', 'decerr', 'radeccov', \
                        'u', 'uerr', 'uesq', 'unobs', 'uerrfl', 'uskyflg', 'ushpflg','uovrflg', \
                        'g', 'gerr', 'gesq', 'gnobs', 'gerrfl', 'gskyflg', 'gshpflg','govrflg', \
                        'r', 'rerr', 'resq', 'rnobs', 'rerrfl', 'rskyflg', 'rshpflg','rovrflg', \
                        'i', 'ierr', 'iesq', 'inobs', 'ierrfl', 'iskyflg', 'ishpflg','iovrflg', \
                        'z', 'zerr', 'zesq', 'znobs', 'zerrfl', 'zskyflg', 'zshpflg','zovrflg', \
                        'y', 'yerr', 'yesq', 'ynobs', 'yerrfl', 'yskyflg', 'yshpflg','yovrflg']

        self.lineHeader=','.join(listHeader)

def go(nsidePow=10):

    """Tester for some healpix functionality. Put nsidePow<=10 to reproduce the NASA HEALPIX frontpage documenation examples."""

    BD = bdbsCat(nsidePow=nsidePow)
    BD.estHealpixSpacing()
    BD.wrapLoadAndID()

def testStream(nMax=-1, nBunch=1000000, splitFiles=True, Verbose=True, infil='TEST.catalog'):

    """Tests stream-based reading of input file. Example calls:

    txt2file.testStream(2000, nBunch=237, splitFiles=True)

    txt2file.testStream(-1, nBunch=100000, splitFiles=True, Verbose=True)"""
    
    BD = bdbsCat(infil, 22)
    BD.setOutputFilename()
    BD.ensureOutdirExists()

    # control the bunch size so we can test our splitter; determine
    # whether output written to screen
    BD.nBunch=nBunch
    BD.Verbose=Verbose

    # are we distributing the output across multiple files?
    BD.splitFiles=splitFiles

    BD.setHeaderLineCSV()

    Debug=False
    if 0 < nMax < 10:
        Debug=True

    BD.processStream(iMax=nMax, DBG=Debug)

