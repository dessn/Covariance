"""Python program to compute C_stat
"""

from optparse import OptionParser

def compute_Cstat(options):
    """Python program to compute C_stat
    """

    import numpy
    import astropy.io.fits as fits
    from astropy.table import Table
    import JLA_library as JLA

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


    # -----------  Read in the data --------------------------

    print 'There are %d SNe in the sample' % (nSNe)

    indices = JLA.reindex_SNe(SNeList['id'], SNe)
    SNe=SNe[indices]

    C_stat=numpy.zeros(9*nSNe*nSNe).reshape(3*nSNe,3*nSNe)

    for i,SN in enumerate(SNe):
        cov=numpy.zeros(9).reshape(3,3)
        cov[0,0]=SN['dmb']**2.
        cov[1,1]=SN['dx1']**2.
        cov[2,2]=SN['dcolor']**2.
        cov[0,1]=SN['cov_m_s']
        cov[0,2]=SN['cov_m_c']
        cov[1,2]=SN['cov_s_c']
        # symmetrise
        cov=cov+cov.T-numpy.diag(cov.diagonal())
        C_stat[i*3:i*3+3,i*3:i*3+3]=cov

    # -----------  Read in the base matrix computed using salt2_stat.cc ------------

    if options.base!=None:
        C_stat+=fits.getdata(options.base)

    date = JLA.get_date()
    fits.writeto('C_stat_%s.fits' % date,C_stat,clobber=True) 

    return
    
if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA parameters")

    PARSER.add_option("-j", "--jla", dest="jla", default=False,
                      action="store_true",
                      help="Restrict this to SNe in the JLA sample")

    PARSER.add_option("-b", "--base", dest="base", default=None,
                      help="Base matrix, at the start, this is the matrix ")

    PARSER.add_option("-s", "--SNlist", dest="SNlist", 
                      help="List of SN")
    
    PARSER.add_option("-l", "--lcfits", dest="lcfits", default="lightCurveFits",
                      help="Key in config file pointing to lightcurve fit parameters")

    PARSER.add_option("-u", "--update", dest="update", default=None,
                      help="A list of the SNe to update")

    (OPTIONS, ARGS) = PARSER.parse_args()

    compute_Cstat(OPTIONS)


