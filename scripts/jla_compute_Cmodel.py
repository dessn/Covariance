"""
Python program to compute C_model
"""

from optparse import OptionParser

# Cmodel accounts for the error that comes from assuming that the SALT2.4 model adequately describes 
# SNe Ia. This is called the light curve model uncertainty. See section 4.4 in B14
# The idea is to simulate SNe data sets according to specific SNe Ia models, train SALT2.4 with a
# training sample derived from that model and then fit a test sample derived from that model.
# They use the models in M14 that lead to the largest Hubble diagram bias. This is the  G10' - C11 model
# There are two approaches to computing the uncertainty here. 
# i) Use the fit to w (and Om) to estimate the offset in the distance modulus, Table 7,
#
# ii) Use Fig 16 directly.

# We use the second approach, as this is the one used in the JLA
# G10'-C11-REAL-REAL => Based on the G10' input model with C11 intrinsic scatter trained 
# This is the bottom left plot in FIg. 16 of Mosher et al. 2014, ApJ, 793, 16


def compute_model(options):

    import numpy
    import astropy.io.fits as fits
    import JLA_library as JLA
    from astropy.table import Table
    from astropy.cosmology import FlatwCDM
    from scipy.interpolate import interp1d


    # -----------  Read in the configuration file ------------
    params=JLA.build_dictionary(options.config)

    # -----------  Read in the SN ordering ------------------------
    SNeList = numpy.genfromtxt(options.SNlist,
                               usecols=(0, 2),
                               dtype='S30,S200',
                               names=['id', 'lc'])
    nSNe = len(SNeList)

    for i, SN in enumerate(SNeList):
        SNeList['id'][i] = SNeList['id'][i].replace('lc-', '').replace('.list', '')

    lcfile = JLA.get_full_path(params[options.lcfits])
    SNe = Table.read(lcfile, format='fits')

    print 'There are %d SNe' % (nSNe)

    # For the JLA SNe
    redshift = SNe['zcmb']
    replace=(redshift < 0)

    # For the non JLA SNe
    redshift[replace]=SNe[replace]['zhel']

    if options.raw:
        # Data from the bottom left hand figure of Mosher et al. 2014.
        # This is option ii) that is descibed above
        offsets=Table.read(JLA.get_full_path(params['modelOffset']),format='ascii.csv')
        Delta_M=interp1d(offsets['z'], offsets['offset'], kind='linear',bounds_error=False,fill_value='extrapolate')(redshift)
    else:
        Om_0=0.303 # JLA value in the wCDM model
        cosmo1 = FlatwCDM(name='SNLS3+WMAP7', H0=70.0, Om0=Om_0, w0=-1.0)
        cosmo2 = FlatwCDM(name='SNLS3+WMAP7', H0=70.0, Om0=Om_0, w0=-1.024)
        Delta_M=5*numpy.log10(cosmo1.luminosity_distance(redshift)/cosmo2.luminosity_distance(redshift))
    
    # Build the covariance matrix. Note that only magnitudes are affected
    Zero=numpy.zeros(nSNe)
    H=numpy.concatenate((Delta_M,Zero,Zero)).reshape(3,nSNe).ravel(order='F')
    C_model=numpy.matrix(H).T * numpy.matrix(H)

    date = JLA.get_date()
    fits.writeto('C_model_%s.fits' % (date),numpy.array(C_model),clobber=True) 

    return None

if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA parameters")

    PARSER.add_option("-s", "--SNlist", dest="SNlist", 
                      help="List of SN")

    PARSER.add_option("-l", "--lcfits", dest="lcfits", default="lightCurveFits",
                      help="Key in config file pointing to lightcurve fit parameters")

    PARSER.add_option("-r", "--raw", dest="raw", default=True,
                      help="Use the raw data points")

    
    (options, args) = PARSER.parse_args()

    compute_model(options)
