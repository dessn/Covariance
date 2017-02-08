"""Python program to compute C_K 
"""


from optparse import OptionParser

# Usage
# JLA_computeCcal.py opitons
# 
# Written specifically for the JLA-like analysis of the DES spectroscopically confirmed sample
# For Y3A1, the calibration consists of the following steps
# Controlled by DES
# 1) The Forward Global Calibration Model (FGCM), see Burke et al.
# 2) Determining the offsets between the standard and synthetic magnitudes of the CALSPEC stars
# 3) The photometric methods used to estimate fluxes for SNe, CALPSEC standards and calibration stars.
# Controlled by STScI
# 4) The calibration of the secondary stanadrds from the primary HST standards
# 5) Definition of the absoluted SEDs of the primary HST standards

# Each step has uncertainties.

def compute_C_K(options):
    import JLA_library as JLA
    import numpy
    import astropy.io.fits as fits

    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)


    # -----------  If required, read in the JLA version of C_Kappa ------------

    nDim=46
    C_K_DES=numpy.zeros(nDim*nDim).reshape(nDim,nDim)
     
    PS1_unc=0.5   # 5 mmag uncertainty in the PanSTARRS calibration
    FGCM_unc=0.5  # 5 mmag uncertainty, see Burke et al.
    nCALSPEC_Observations=100
    AB_unc=0.5 / numpy.sqrt (nCALSPEC_Observations)

    if options.base:
        # Read in the JLA matrix and extract the appropriate rows and columns
        # The matrix diagonal contains uncertainties in the ZPs first,
        # and uncertainties in the filter curves second
        # The order is specified in salt2_calib_variations_all/saltModels.list
        # CfA3 and CfA4 rows 10 to 19. We write these to rows 1 to 10
        # CSP  rows 20 to 25. We write these to rows 11 to 15
        C_K_JLA=fits.getdata(JLA.get_full_path(params['C_kappa_JLA']))

        # Extract the relevant columns and rows
        # ZPs first
        size=C_K_JLA.shape[0]
        sel=numpy.zeros(size,bool)
        sel[9:19]=True
        sel[19:25]=True
        sel2d= numpy.matrix(sel).T * numpy.matrix(sel)
        C_K_DES[0:16,0:16]=C_K_JLA[sel2d].reshape(16,16)
        
        # Filter curves second
        sel=numpy.zeros(size,bool)
        sel[9+size/2:19+size/2]=True
        sel[19+size/2:25+size/2]=True
        sel2d= numpy.matrix(sel).T * numpy.matrix(sel)
        C_K_DES[23:39,23:39]=C_K_JLA[sel2d].reshape(16,16)

    # Compute the terms in DES, this includes the cross terms as well
    # We first compute them separately, then add them to the matrix
    # Read in the table listing the uncertainties in the ZPs and effective wavelengths


    filterUncertainties=numpy.genfromtxt(JLA.get_full_path(params['filterUncertainties']),comments='#',usecols=(0,1,2,3,5,6),dtype='S30,f8,f8,f8,f8,f8',names=['filter','zp','wavelength','central','relative_ZP','filterFactor'])

    nFilters=len(filterUncertainties)
    C_K_new=numpy.zeros(nFilters*nFilters*4).reshape(nFilters*2,nFilters*2)
    #1) and #2) The variance from the uncertainties in the ZP and the central wavelengths

    for i,filt in enumerate(filterUncertainties):
        C_K_new[i,i]=(filt['zp']/1000.)**2. - (PS1_unc/1000.)**2. + (FGCM_unc/1000.)**2. + (AB_unc/1000.)**2.
        C_K_new[i+nFilters,i+nFilters]=(filt['filterFactor']*filt['wavelength'])**2.
        

    #4) B14 3.4.1 The uncertainty in the colour of the WD system 0.5% from 3,000-10,000
    # The uncertainty is computed with respect to the Bessell B filter.
    # We are not concerned about the absoulte uncertainty as this is degenegate with H_0 and the absolute magnitude of the SNe.

    slope=0.005
    waveStart=300.
    waveEnd=1000.
    B_central=436.0

    for i,filt1 in enumerate(filterUncertainties):
        for j,filt2 in enumerate(filterUncertainties):
            if i>=j:
#                C_K_new[i,j]+=(slope / (waveEnd - waveStart) * (filt1['central']-filt2['central']))**2.
                C_K_new[i,j]+=(slope / (waveEnd - waveStart) * (filt1['central']-B_central)) * \
                              (slope / (waveEnd - waveStart) * (filt2['central']-B_central))

    # We need to include the cross terms between DES and the nearby SNe

    #3) B14 3.4.1 The uncertainty associated to the measurement of the Secondary CALSPEC standards
    # The uncerteinty is assumed to be uncorrelated between filters
    # It only affects the diagonal terms of the ZPs
    # It is estmated from repeat STIS measurements of the standard AGK+81D266  Bohlin et al. 2000 AJ 120, 437


    nObs_C26202=1          # It's been observed at least once
    unc_transfer=0.003     # See Bohlin et al. 2000 AJ 120, 437 and Bohlin 1999 ISR 99-07
    for i,filt1 in enumerate(filterUncertainties):
        C_K_new[i,i]+=unc_transfer**2. / nObs_C26202

    C_K_new=C_K_new+C_K_new.T-numpy.diag(C_K_new.diagonal())
    
    # Update C_K. We do not update the terms that come from JLA
    # For the Bc filter of CfA, and the V1 and V2 filters of CSP, we asumme that they have 
    # the same sized systematic uncertainteies as B filter of CfA and V1 and V2 filters of CSP
    sel=numpy.zeros(nDim,bool)
    sel[0:16]=True
    sel[23:39]=True
    sel2d= numpy.matrix(sel).T * numpy.matrix(sel)
    C_K_new[sel2d]=0.0
    C_K_DES+=C_K_new

    # Write out the results
    date=JLA.get_date()
    hdu = fits.PrimaryHDU(C_K_DES)
    hdu.writeto("%s_%s.fits" % (options.output,date),clobber=True)

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-b", "--base", dest="base", default=True,
                      action="store_true",
                      help="Use the JLA matrix as the base")

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA files")

    parser.add_option("-o", "--output", dest="output", default="DES_C_K",
                      help="Output file")

    (options, args) = parser.parse_args()


    compute_C_K(options)
