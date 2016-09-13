"""Python program to compute C_dust for a new SNe
"""

from optparse import OptionParser

def compute_dust(options):
    """Python program to compute C_dust
    """

    import numpy
    import astropy.io.fits as fits
    import os
    import JLA_library as JLA

    # ---------- Read in the SNe list -------------------------

    SNelist = numpy.genfromtxt(options.SNlist,
                               usecols=(0, 2),
                               dtype='S30,S100',
                               names=['id', 'lc'])

    for i, SN in enumerate(SNelist):
        SNelist['id'][i] = SNelist['id'][i].replace('lc-','').replace('.list','')

    # -----------   The lightcurve fitting -------------------

    # Compute the offset between the nominal value of the extinciton 
    # and the adjusted value
    offset = 0.1

    j = []

    for SN in SNelist:
        inputFile = SN['lc']
        print 'Fitting %s' % (SN['id'])
        dm, dx1, dc = JLA.compute_extinction_offset(SN['id'], inputFile, offset, options.workArea)
        j.extend([dm, dx1, dc])
    
    cdust = numpy.matrix(j).T * numpy.matrix(j) * 4.0

    date = JLA.get_date()

    fits.writeto('C_dust_%s.fits' % date, cdust, clobber=True) 

    return

if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", default="JLA.config",
                  help="Parameter file containting the location of various JLA parameters")

    PARSER.add_option("-s", "--SNlist", dest="SNlist", 
                  help="List of SN")
        
    PARSER.add_option("-w", "--workArea", dest="workArea", default="../workArea",
                      help="Work area that stores the light curve fits")

    (OPTIONS, ARGS) = PARSER.parse_args()

    compute_dust(OPTIONS)

