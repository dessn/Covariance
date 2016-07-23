"""Python program to compute the magnitude of the chosen CALSPEC standard
"""

from optparse import OptionParser

def compute_ZP(options):

    import JLA_library as JLA
    import numpy as np

    params=JLA.build_dictionary(options.config)

    # Read in the standard star

    standard=JLA.spectrum(JLA.get_full_path(params['magSys'])+options.standard)

    # Read in the filter

    filt=JLA.filterCurve(JLA.get_full_path(params['filterDir'])+options.filter)

    # Compute the ZP
    if options.system=='AB':
        print '%s in %s %s %5.3f' % (options.standard,options.filter,options.system,filt.AB(standard))
    else:
        pass
#        print '%s in %s %s %5.3f' % (options.standard,options.filter,options.system,filt.Vega(standard))


    return


if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="DES.config",
                      help="Parameter file containting the location of various JLA files")

    parser.add_option("-s", "--standard", dest="standard", default=None,
                      help="Standard Star Name")

    parser.add_option("-f", "--filter", dest="filter", default=None,
                      help="Filter")

    parser.add_option("-S", "--system", dest="system", default="AB",
                      help="System (AB/Vega)")
    
    (options, args) = parser.parse_args()

    compute_ZP(options)
