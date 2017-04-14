"""
Python program to compute C_K
"""

from optparse import OptionParser

# Usage
# JLA_computeCcal.py -c config_file
# The config file contains the locations of files needed by the script
#
# Written specifically for the JLA-like analysis of the DES spectroscopically
# confirmed sample
# For Y3A1, the calibration consists of the following steps
# Those controlled by DES
# 1) The Forward Global Calibration Model (FGCM), see Burke et al.
# 2) The offsets between the FGCM and the AB system.
# 3) The methods used to estimate fluxes for SNe, CALPSEC standards
# and calibration stars. This is set to zero as the methods are the same.
# Those controlled by STScI
# 4) The calibration of the secondary standards from the primary HST standards
# 5) Definition of the absoluted SEDs of the primary HST standards

# Each step has uncertainties.


def compute_C_K(options):
    import JLA_library as JLA
    import numpy
    import astropy.io.fits as fits

    # -----------  Read in the configuration file ------------

    params = JLA.build_dictionary(options.config)

    # -----------  We read in the JLA version of C_Kappa ------------

    nDim = 46  # The number of elements in the DES C_Kappa matrix
    C_K_DES = numpy.zeros(nDim * nDim).reshape(nDim, nDim)

    PS1_unc = 0.005   # 5 mmag uncertainty in the PanSTARRS calibration
    FGCM_unc =0.005  # 5 mmag uncertainty, see Burke et al.
    nC26202_Observations = 100  # Needs revision
    AB_unc = FGCM_unc / numpy.sqrt(nC26202_Observations)

    if options.base:
        # Read in the JLA matrix and extract the appropriate rows and columns
        # The matrix is structured in blocks with ZPs first,
        # and uncertainties in the filter curves second
        # The order is specified in salt2_calib_variations_all/saltModels.list
        # CfA3 and CfA4 are in rows 10 to 19 and 46 to 56 (starting at row 1) 
        # We write these to rows 1 to 10 and 24 to 33
        # CSP are in rows 20 to 25 and 57 to 62.
        # We write these to rows 11 to 16 and 34 to 39
        
        C_K_JLA = fits.getdata(JLA.get_full_path(params['C_kappa_JLA']))

        # Extract the relevant columns and rows
        # ZPs first
        size = C_K_JLA.shape[0]
        sel = numpy.zeros(size, bool)
        sel[9:19] = True
        sel[19:25] = True
        sel2d = numpy.matrix(sel).T * numpy.matrix(sel)
        C_K_DES[0:16, 0:16] = C_K_JLA[sel2d].reshape(16, 16)

        # Filter curves second
        sel = numpy.zeros(size, bool)
        sel[9+size/2:19+size/2] = True
        sel[19+size/2:25+size/2] = True
        sel2d = numpy.matrix(sel).T * numpy.matrix(sel)
        C_K_DES[23:39, 23:39] = C_K_JLA[sel2d].reshape(16, 16)

    # Read in the table listing the uncertainties in the ZPs and
    # effective wavelengths

    filterUncertainties = numpy.genfromtxt(JLA.get_full_path(params['filterUncertainties']),
                comments='#',usecols=(0,1,2,3,5,6), dtype='S30,f8,f8,f8,f8,f8',
                names=['filter', 'zp', 'wavelength', 'central', 'relative_ZP', 'filterFactor'])

    # For the Bc filter of CfA, and the V1 and V2 filters of CSP,
    # we asumme that they have the same sized systematic uncertainteies as
    # B filter of CfA and V1 and V2 filters of CSP
    # We could either copy these terms across or recompute them.
    # We choose to recompute them

    # Compute the terms in DES, this includes the cross terms
    # We first compute them separately, then add them to the matrix

    nFilters = len(filterUncertainties)
    C_K_new = numpy.zeros(nFilters*nFilters*4).reshape(nFilters*2, nFilters*2)

    # 1 and 2) The contribution from the uncertainties in the ZPs and the
    #   central wavelengths of the filter curves
    #   This uncertainty in the ZP has seeral components
    #   a) The SN field-to-field variation between DES and PS1
    #   b) PS1 is not perfect, hence we reduce the uncertainty
    #   c) The uncertainty in the measurement of the transfer to the AB system
    #      using the observations of C26202
    # Not included yet is the systematic uncertainty that comes from
    # chromatic correction
    # We do not include FGCM_unc as this is a statistical error that
    # should be included in the lightcurve uncetainties

    # We compute all terms but insert the JLA computed terms later

    for i, filt in enumerate(filterUncertainties):
        C_K_new[i, i] = (filt['zp'] / 1000.)**2. - PS1_unc**2. + AB_unc**2.
        C_K_new[i+nFilters, i+nFilters] = (filt['wavelength'])**2.

    # 3) Set to zero for DES

    # 5) B14 3.4.1 The uncertainty in the colour of the WD system 0.5%
    # from 3,000-10,000
    # The uncertainty is computed with respect to the Bessell B filter.
    # The Bessell B filter is the filter we use in computing the dist. modulus
    # The absolute uncertainty at the rest frame wavelengt of the B band
    # is not important here, as this is absorbed into the
    # combination of the absolute B band magnitude of SNe Ia and
    # the Hubble constant.

    slope = 0.005
    waveStart = 300.
    waveEnd = 1000.
    B_central = 436.0

    for i, filt1 in enumerate(filterUncertainties):
        for j, filt2 in enumerate(filterUncertainties):
            if i >= j:
                C_K_new[i, j] += (slope / (waveEnd - waveStart) * (filt1['central']-B_central)) * \
                                 (slope / (waveEnd - waveStart) * (filt2['central']-B_central))

    # 4) B14 3.4.1 The uncertainty associated to the measurement of
    # the Secondary CALSPEC standards
    # The uncerteinty is assumed to be uncorrelated between filters
    # It only affects the diagonal terms of the ZPs
    # It is estmated from repeat STIS measurements of the standard
    # AGK+81D266  Bohlin et al. 2000 AJ 120, 437 and Bohlin 1999 ISR 99-07

    nObs_C26202 = 1          # It's been observed once
    unc_transfer = 0.003

    for i, filt1 in enumerate(filterUncertainties):
        C_K_new[i, i] += unc_transfer**2. / nObs_C26202

    C_K_new = C_K_new+C_K_new.T-numpy.diag(C_K_new.diagonal())

    # Update C_K. We do not update the terms that come from JLA

    sel = numpy.zeros(nDim, bool)
    sel[0:16] = True
    sel[23:39] = True
    sel2d = numpy.matrix(sel).T * numpy.matrix(sel)
    C_K_new[sel2d] = 0.0
    C_K_DES += C_K_new

    # Write out the results
    date = JLA.get_date()
    hdu = fits.PrimaryHDU(C_K_DES)
    hdu.writeto("%s_%s.fits" % (options.output, date), clobber=True)

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-b", "--base", dest="base", default=True,
                      action="store_true",
                      help="Use the JLA matrix as the base")

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the \
                          location of various JLA files")

    parser.add_option("-o", "--output", dest="output", default="DES_C_K",
                      help="Output file")

    (options, args) = parser.parse_args()

    compute_C_K(options)
