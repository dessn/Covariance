"""Python that estimates the change in the effective wavelength of the DES filters across the DES focal plane
"""

from optparse import OptionParser

# Usage
# JLA_computeCcal.py opitons
# 
# Written specifically for the JLA-like analysis of the DES spectroscopically confirmed sample
# Includes the relevant columns from the 

def compute_filterTransVar(options):
    import JLA_library as JLA
    import numpy
    from astropy.table import Table

    eff=[]
    # Read in the filter curves
    filt=Table.read(JLA.get_full_path(options.filter),format='ascii.csv')
    for col in filt.colnames:
        if 'ccd' in col and 'amp' not in col:
            f=JLA.filterCurve(filt['wavelength'],filt[col])
            eff.append(f.eff())

    print 'Examined the transmission curves for %d CCDs' % (len(eff))
    print 'Mean effective wavelength is %6.1f' % (numpy.mean(eff))
    print 'RMS effective wavelength %6.1f' % (numpy.std(eff))

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-f", "--filter", dest="filter", default=None,
                      help="Input Filter")

    (options, args) = parser.parse_args()

    compute_filterTransVar(options)
