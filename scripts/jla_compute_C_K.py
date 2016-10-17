"""Python program to compute C_K 
"""

from optparse import OptionParser

# Usage
# JLA_computeCcal.py opitons
# 
# Written specifically for the JLA-like analysis of the DES spectroscopically confirmed sample
# Includes the relevant columns from the 

def compute_C_K(options):
    import JLA_library as JLA
    import numpy
    import astropy.io.fits as fits

    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)


    # -----------  If required, read in the JLA version of C_Kappa ------------

    nDim=36
    C_K_DES=numpy.zeros(nDim*nDim).reshape(nDim,nDim)
     

    if options.base:
        # Read in the JLA matrix and extract the apprpriate rows and columns
        # The matrix diagonal contains uncertainties in the ZPs first,
        # and uncertainties in the filter curves second
        # The order is specified in salt2_calib_variations_all/saltModels.list
        # CfA3 and CfA4 rows 10 to 14 We write these to rows 1 
        # CSP  rows 20 to 25
        C_K_JLA=fits.getdata(JLA.get_full_path(params['C_kappa_JLA']))

        # Extract the relevant columns and rows
        # ZPs first
        size=C_K_JLA.shape[0]
        sel=numpy.zeros(size,bool)
        sel[9:14]=True
        sel[19:25]=True
        sel2d= numpy.matrix(sel).T * numpy.matrix(sel)
        C_K_DES[0:11,0:11]=C_K_JLA[sel2d].reshape(11,11)
        # We then add Bc, V1 and V2 - to do
        
        # Filter curves second
        sel=numpy.zeros(size,bool)
        sel[9+size/2:14+size/2]=True
        sel[19+size/2:25+size/2]=True
        sel2d= numpy.matrix(sel).T * numpy.matrix(sel)
        C_K_DES[18:29,18:29]=C_K_JLA[sel2d].reshape(11,11)

    # Compute the terms in DES, this includes the cross terms as well
    # We first compute them separately, then add them to the matrix
    # Read in the table listing the uncertainties in the ZPs and effective wavelengths

    filterUncertainties=numpy.genfromtxt(JLA.get_full_path(params['filterUncertainties']),comments='#',usecols=(0,1,2,3,5),dtype='S30,f8,f8,f8,f8',names=['filter','zp','wavelength','central','relative_ZP'])

    nFilters=len(filterUncertainties)
    C_K_new=numpy.zeros(nFilters*nFilters*4).reshape(nFilters*2,nFilters*2)

    # The variance from the uncertainties in the ZP and the central wavelengths

    for i,filt in enumerate(filterUncertainties):
        C_K_new[i,i]=(filt['zp']/1000.)**2.
        C_K_new[i+nFilters,i+nFilters]=filt['wavelength']**2.
        

    # B14 3.4.1 The uncertainty in the colour of the WD system 0.5% from 3,000-10,000

    slope=0.005
    waveStart=300.
    waveEnd=1000.

    for i,filt1 in enumerate(filterUncertainties):
        for j,filt2 in enumerate(filterUncertainties):
            if i>=j:
                C_K_new[i,j]+=(slope / (waveEnd - waveStart) * (filt1['central']-filt2['central']))**2.


    # We need to include the cross terms between DES and the nearby SNe

    # B14 3.4.1 The uncertainty associated to the measurement of the Secondary CALSPEC standards
    # This only affects the diagonal terms of the ZPs
                
    # Uncertainties specific to the calibration of DES
    # See the steps that are involved in the GRC model
    # The most important steps are        
    # i) The global relative calibration, which is applied to single exposures
    # The SN fields used the individual standard star fields a anchors.
    # This will be different for Y3A1

    for i,filt in enumerate(filterUncertainties):
        C_K_new[i,i]+=(filt['relative_ZP']/1000.)**2.

    # ii) The SLR adjustment
    # It is not clear if this is applied to the single exposures. 
    # We assume that it is not applied to single exposures

    C_K_new=C_K_new+C_K_new.T-numpy.diag(C_K_new.diagonal())
    
    # Update C_K. We do not update the terms that come from JLA
    # For the Bc filter of CfA, and the V1 and V2 filters of CSP, we asumme that they have the same sized systematic uncertainteies
    # as B filter of CfA and V1 and V2 filters of CSP
    sel=numpy.zeros(nDim,bool)
    sel[0:11]=True
    sel[18:29]=True
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
