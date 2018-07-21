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
# See 

def compute_C_K(options):
    import JLA_library as JLA
    import jla_FGCM as FGCM
    import numpy
    import astropy.io.fits as fits

    # -----------  Read in the configuration file ------------

    params = JLA.build_dictionary(options.config)

    # -----------  We read in the JLA version of C_Kappa ------------

    nDim = 52  # The number of elements in the DES C_Kappa matrix
    C_K_DES = numpy.zeros(nDim * nDim).reshape(nDim, nDim)

    SMP_ZP = 0.001                                 # The accuracy of the SMP ZPs
    # 6.6 mmag RMS scatter between FGCM and GAIA
    # The first factor of sqrt(2) comes from asuuming that GAIA and DES contrinute eually to the error
    # The second factor of sqrt(2) arises because we are taking the difference between two points, one where the standard is, and another where the SN is.
    FGCM_GAIA = 0.0066 / numpy.sqrt(2.) * numpy.sqrt(2.)    
    nC26202_Observations = {'DES_g':133,'DES_r':21,'DES_i':27,'DES_z':78}                      # Number of times C26202 has been observed
    FGCM_unc = 0.005                               # The RMS scatter in FGCM standard magnitudes
    chromatic_differential = 0.0                   # Set to zero for now

    if options.base:
        # Read in the JLA matrix and extract the appropriate rows and columns
        # The matrix is structured in blocks with ZPs first,
        # and uncertainties in the filter curves second
        # The order is specified in salt2_calib_variations_all/saltModels.list
        # CfA3 and CfA4 are in rows 10 to 19 and 46 to 56 (starting at row 1) 
        # We write these to rows 1 to 10 and 27 to 36
        # CSP are in rows 20 to 25 and 57 to 62.
        # We write these to rows 11 to 16 and 37 to 42
        
        C_K_JLA = fits.getdata(JLA.get_full_path(params['C_kappa_JLA']))

        # Extract the relevant columns and rows
        # ZPs first
        # Since the indices for CfA4, CfA4, and CSP are consecutive, we do this all at once
        size = C_K_JLA.shape[0]
        C_K_DES[0:16, 0:16] = C_K_JLA[9:25,9:25]

        # Filter curves second
        C_K_DES[27:42, 27:42] = C_K_JLA[9+size/2:25+size/2,9+size/2:25+size/2]

        # Cross terms. Not needed, as they are zero
        # C_K_DES[0:16, 23:39] = C_K_JLA[9:25,9+size/2:25+size/2]
        # C_K_DES[23:39, 0:16] = C_K_JLA[9+size/2:25+size/2,9:25]

    # Read in the table listing the uncertainties in the ZPs and
    # effective wavelengths

    filterUncertainties = numpy.genfromtxt(JLA.get_full_path(params['filterUncertainties']),
                comments='#',usecols=(0,1,2,3,4), dtype='S30,f8,f8,f8,f8',
                names=['filter', 'zp', 'zp_off', 'wavelength', 'central'])

    # For the Bc filter of CfA, and the V1 and V2 filters of CSP,
    # we asumme that they have the same sized systematic uncertainteies as
    # B filter of CfA and V1 and V2 filters of CSP
    # We could either copy these terms across or recompute them.
    # We choose to recompute them 


    # Compute the terms in DES, this includes the cross terms
    # We first compute them separately, then add them to the matrix

    nFilters = len(filterUncertainties)
    C_K_new = numpy.zeros(nFilters*nFilters*4).reshape(nFilters*2, nFilters*2)

    # 1) DES controlled uncertainties  
    #   This uncertainty in the ZP has seeral components
    #   a) The uncertainty in the differential chromatic correction (set to zero for now)
    #   Note that this error is 100% correlated to the component of b) that comes from the filter curve
    #   b) The uncertainty in the measurement of the transfer to the AB system
    #      using the observations of C26202
    #   c) The SN field-to-field variation between DES and GAIA

    for i, filt in enumerate(filterUncertainties):
        if 'DES' in filt['filter']:
            error_I0,error_chromatic,error_AB=FGCM.prop_unc(params,filt)
            #print numpy.sqrt((error_AB)**2. + (FGCM_unc)**2. / nC26202_Observations[filt['filter']])
            C_K_new[i, i] = FGCM_GAIA**2. + (error_AB)**2.+(FGCM_unc)**2. / nC26202_Observations[filt['filter']] + SMP_ZP**2.
            print '%s %5.4f' % (filt['filter'],numpy.sqrt(C_K_new[i, i]))
            C_K_new[i, i+nFilters] = (error_AB) * filt['wavelength']
            C_K_new[i+nFilters, i] = (error_AB) * filt['wavelength']
            C_K_new[i+nFilters, i+nFilters] = (filt['wavelength'])**2.
        else:
            C_K_new[i, i] = (filt['zp'] / 1000.)**2. + (filt['zp_off'] / 3. / 1000.)**2.


    # 2a) B14 3.4.1 The uncertainty associated to the measurement of
    # the Secondary CALSPEC standards
    # The uncerteinty is assumed to be uncorrelated between filters
    # It only affects the diagonal terms of the ZPs
    # It is estmated from repeat STIS measurements of the standard
    # AGK+81D266  Bohlin et al. 2000 AJ 120, 437 and Bohlin 1999 ISR 99-07

    nObs_C26202 = 1          # It's been observed once
    unc_transfer = 0.003     # 0.3% uncertainty

    for i, filt1 in enumerate(filterUncertainties):
        C_K_new[i, i] += unc_transfer**2. / nObs_C26202


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
                C_K_new[i, j] += (slope / (waveEnd - waveStart) * (filt1['central']-central)) * \
                                 (slope / (waveEnd - waveStart) * (filt2['central']-central))
                
    C_K_new = C_K_new+C_K_new.T-numpy.diag(C_K_new.diagonal())

    if options.base:
        # We do not update 
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

    parser.add_option("-b", "--base", dest="base", default=False,
                      action="store_true",
                      help="Use the JLA matrix as the base")

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the \
                          location of various JLA files")

    parser.add_option("-o", "--output", dest="output", default="DES_C_K",
                      help="Output file")

    (options, args) = parser.parse_args()

    compute_C_K(options)
