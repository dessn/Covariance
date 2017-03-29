"""Takes the filter curves that Ting Li provided and saves them in a format that can be used by SALT2
"""

from optparse import OptionParser
from astropy.table import Table
import JLA_library as JLA

bounds={'g':{"lower":3800,"upper":5580},
        'r':{"lower":5290,"upper":7320},
        'i':{"lower":6710,"upper":8910},
        'z':{"lower":8050,"upper":10385}}

def make_FilterCurves(options):

    filterCurve=Table.read(options.input,format='fits')

    f=open("des_y3a1_std_%s.dat" % (options.filterName),'w')
    f.write("# Written on %s\n" % (JLA.get_date()))
    f.write("# Derived from %s\n" % (options.input))
    f.write("# Wavelength (Angstroms) Transmission\n")
    selection=(filterCurve["lambda"] >  bounds[options.filterName]["lower"]) & (filterCurve["lambda"] <  bounds[options.filterName]["upper"])
    for line in filterCurve[selection]:
        f.write("%5.1f %7.5f\n" % (line["lambda"],line[options.filterName]))

    f.close

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="input", default=None,
                      help="input filter filename")

    parser.add_option("-f", "--filterName", dest="filterName", default=None,
                      help="filter name")

    (options, args) = parser.parse_args()


    make_FilterCurves(options)
