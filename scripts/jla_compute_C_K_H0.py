"""
Python program to compute C_K
"""

from optparse import OptionParser

# Usage
# JLA_computeCcal.py -c config_file
# The config file contains the locations of files needed by the script
# Written specifically for the JLA-like analysis of Bonnie's H0 sample


def compute_C_K(options):
    import JLA_library as JLA
    import numpy
    import astropy.io.fits as fits

    # -----------  Read in the configuration file ------------

    params = JLA.build_dictionary(options.config)

    # -----------  We read in the JLA version of C_Kappa ------------

    if options.base:
        # CfA1 and CfA2 not treated separately and we use the JLA uncertainties
        nDim = 42
    else:
        # CfA1 and CfA2 treated separately, and we use the Pantheon uncertainties
        nDim = 58
        
    C_K_H0 = numpy.zeros(nDim * nDim).reshape(nDim, nDim)

    if options.base:
        # Read in the JLA matrix and extract the appropriate rows and columns
        # The matrix is structured in blocks with ZPs first,
        # and uncertainties in the filter curves second
        # The order is specified in salt2_calib_variations_all/saltModels.list
        # Standard, Landolt photometry is in rows 5 to 9 and rows 42 to 46
        # Keplercam is in rows 10 to 14 and 47 to 51
        # 4 Shooter is in rows 15 to 19 and 52 5o 56
        # CSP is in rows 20 to 25 and 56 to 62
        
        C_K_JLA = fits.getdata(JLA.get_full_path(params['C_kappa_JLA']))

        # Extract the relevant columns and rows
        # ZPs first
        # Since the indices are consecutive, we do this all at once
        size = C_K_JLA.shape[0]
        C_K_H0[0:21, 0:21] = C_K_JLA[4:25,4:25]

        # Filter curves second
        C_K_H0[21:42, 21:42] = C_K_JLA[4+size/2:25+size/2,4+size/2:25+size/2]
    else:
        filterUncertainties = numpy.genfromtxt(JLA.get_full_path(params['filterUncertainties']),
                comments='#',usecols=(0,1,2,3,4), dtype='S30,f8,f8,f8,f8',
                names=['filter', 'zp', 'zp_off', 'wavelength', 'central'])

        # 1) ZP and filter uncertainty
        # We add a third of the offset found in Scolnic et al.
        for i, filt in enumerate(filterUncertainties):
            C_K_H0[i, i] = (filt['zp'] / 1000.)**2. + (filt['zp_off'] / 3. / 1000.)**2.
            C_K_H0[i+29, i+29] = (filt['wavelength'])**2.


        # 2a) B14 3.4.1 The uncertainty associated to the measurement of
        # the Secondary CALSPEC standards
        # The uncerteinty is assumed to be uncorrelated between filters
        # It only affects the diagonal terms of the ZPs
        # It is estmated from repeat STIS measurements of the standard
        # AGK+81D266  Bohlin et al. 2000 AJ 120, 437 and Bohlin 1999 ISR 99-07

        # This is the most pessimistic option. We assume that only one standard was observed
        nObs = 1                 # It's been observed once
        unc_transfer = 0.003     # 0.3% uncertainty

        for i, filt1 in enumerate(filterUncertainties):
            C_K_H0[i, i] += unc_transfer**2. / nObs


        # 2b) B14 3.4.1 The uncertainty in the colour of the WD system 0.5%
        # from 3,000-10,000
        # The uncertainty is computed with respect to the Bessell B filter.
        # The Bessell B filter is the filter we use in computing the dist. modulus
        # The absolute uncertainty at the rest frame wavelengt of the B band
        # is not important here, as this is absorbed into the
        # combination of the absolute B band magnitude of SNe Ia and
        # the Hubble constant.

        slope = 0.005
        waveStart = 300
        waveEnd = 1000.
        # central = 436.0 # Corresponds to B filter
        central = 555.6   # Used in the Pantheon sample

        # Note that 2.5 * log_10 (1+x) ~ x for |x| << 1
        for i, filt1 in enumerate(filterUncertainties):
            for j, filt2 in enumerate(filterUncertainties):
                if i >= j:
                    C_K_H0[i, j] += (slope / (waveEnd - waveStart) * (filt1['central']-central)) * \
                        (slope / (waveEnd - waveStart) * (filt2['central']-central))
                

        C_K_H0 = C_K_H0+C_K_H0.T-numpy.diag(C_K_H0.diagonal())


    # Write out the results
    date = JLA.get_date()
    hdu = fits.PrimaryHDU(C_K_H0)
    hdu.writeto("%s_%s.fits" % (options.output, date), clobber=True)

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-b", "--base", dest="base", default=False,
                      action="store_true",
                      help="Use the JLA matrix as the base")

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the \
                          location of various JLA files")

    parser.add_option("-o", "--output", dest="output", default="H0_C_K",
                      help="Output file")

    (options, args) = parser.parse_args()

    compute_C_K(options)
