"""
Python program to estimate the unceetainty in the ZP between two CALPSEC standards
"""

from optparse import OptionParser


def compute_CALPSEC_unc(options):
    import JLA_library as JLA
    import jla_FGCM as FGCM
    import numpy
    import astropy.io.fits as fits

    # -----------  Read in the configuration file ------------

    params = JLA.build_dictionary(options.config)

    filterUncertainties = numpy.genfromtxt(JLA.get_full_path(params['filterUncertainties']),
                comments='#',usecols=(0,1,2,3,4), dtype='S30,f8,f8,f8,f8',
                names=['filter', 'zp', 'zp_off', 'wavelength', 'central'])

    # Compute the magnitude offset for each filter and each spectrum
    for filt in filterUncertainties:
        if 'DES' in filt['filter']:
            error_I0,error_chromatic,error_AB_1=FGCM.prop_unc(params,filt,options.standard1)
            error_I0,error_chromatic,error_AB_2=FGCM.prop_unc(params,filt,options.standard2)
            print '%s: %5.4f' % (filt['filter'],error_AB_1 - error_AB_2)

    return

if __name__ == '__main__':

    parser = OptionParser()


    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the \
                          location of various JLA files")

    parser.add_option("-1", "--standard1", dest="standard1", default=None,
                      help="Path to 1st standard")

    parser.add_option("-2", "--standard2", dest="standard2", default=None,
                      help="Path to 2nd standard")

    (options, args) = parser.parse_args()

    compute_CALPSEC_unc(options)
