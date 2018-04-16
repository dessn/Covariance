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
                               dtype='S30,S110',
                               names=['id', 'lc'])

    for i, SN in enumerate(SNelist):
        SNelist['id'][i] = SNelist['id'][i].replace('lc-','').replace('.list','')

    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)
    try:
        salt_path = JLA.get_full_path(params['defsaltModel'])
    except KeyError:
        salt_path = ''
        
    # -----------   The lightcurve fitting -------------------

    # Compute the offset between the nominal value of the extinciton 
    # and the adjusted value
    # We first compute the difference in light curve fit parameters for E(B-V) * (1+offset)
    offset = 0.1

    j = []

    for SN in SNelist:
        inputFile = SN['lc']
        print 'Fitting %s ' % (SN['lc'])
        workArea = JLA.get_full_path(options.workArea)
        dm, dx1, dc = JLA.compute_extinction_offset(SN['id'], inputFile, offset, workArea, salt_path)
        j.extend([dm, dx1, dc])
    
    # But we want to compute the impact of an offset that is twice as large, hence the factor of 4 in the expression
    # 2017/10/13
    # But we want to compute the impact of an offset that is half as large, hence the factor of 4 in the denominator
    # cdust = numpy.matrix(j).T * numpy.matrix(j) * 4.0
    cdust = numpy.matrix(j).T * numpy.matrix(j) / 4.0

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

