"""Pythom program to compute the systematic unsertainty related to
the contamimation from Ibc SNe"""

# 20161121
# Revised to accept input for the JLA++host sample

# 20170419
# Revised for DES

# 20170504
# remove -1 in digitized (bug) BZ


from optparse import OptionParser

def compute_nonIa(options):
    """Pythom program to compute the systematic unsertainty related to
    the contamimation from Ibc SNe"""

    import numpy
    import astropy.io.fits as fits
    from astropy.table import Table, MaskedColumn, vstack
    import JLA_library as JLA

    # The program computes the covaraince for the spectroscopically confirmed SNe Ia only
    # The prgram assumes that the JLA SNe are first in any list
    # Taken from C11

    # Inputs are the rates of SNe Ia and Ibc, the most likely contaminant

    # Ia rate - Perett et al.
    # SN Ibc rate - proportional to the star formation rate - Hopkins and Beacom
    # SN Ib luminosity distribution. Li et al + bright SN Ibc Richardson

    # The bright Ibc population
    # d_bc = 0.25     # The offset in magnitude between the Ia and bright Ibc
    # s_bc = 0.25     # The magnitude spread
    # f_bright = 0.25 # The fraction of Ibc SN that are bright

    # Simulate the characteristics of the SNLS survey
    # Apply outlier rejection
    # All SNe that pass the cuts are included in the sample

    # One then has a mixture of SNe Ia and SNe Ibc
    # and the average magnitude at each redshift is biased. This
    # is called the raw bias. One multiplies the raw bias by the fraction of
    # objects classified as SNe Ia*

    # The results are presented in 7 redshift bins defined in table 14 of C11
    # We use these results to generate the matrix.
    # Only the SNLS SNe in the JLA sample are considered.
    # For the photometrically selected sample and other surveys, this will probably be different
    # JLA compute this for the SNLS sample only

    # We assume that the redshift in this table refers to the left hand edge of each bin

    # -----------  Read in the configuration file ------------

    params = JLA.build_dictionary(options.config)
    

    data=numpy.genfromtxt(JLA.get_full_path(params['classification']),comments="#",usecols=(0,1,2),dtype=['float','float','float'],names=['redshift','raw_bias','fraction'])
    z_bin=data['redshift']
    raw_bias=data['raw_bias']
    f_star=data['fraction']
    
    # The covaraiance between SNe Ia in the same redshift bin is fully correlated
    # Otherwise, it is uncorrelated

    # -----------  Read in the configuration file ------------

    params = JLA.build_dictionary(options.config)

    SNeList = numpy.genfromtxt(options.SNlist,
                               usecols=(0, 2),
                               dtype='S30,S200',
                               names=['id', 'lc'])

    for i, SN in enumerate(SNeList):
        SNeList['id'][i] = SNeList['id'][i].replace('lc-', '').replace('.list', '')

    lcfile = JLA.get_full_path(params[options.lcfits])
    SNe = Table.read(lcfile, format='fits')

    # Add a bin column and a column that specifies if the covariance needs to be computed
    SNe['bin'] = 0
    SNe['eval'] = False

    # make order of data (in SNe) match SNeList
    
    indices = JLA.reindex_SNe(SNeList['id'], SNe)
    SNe = SNe[indices]

    nSNe = len(SNe)
    # Identify the SNLS SNe in the JLA sample
    # We use the source and the name to decide if we want to add corrections for non-Ia contamination
    # Identify the DESS SNe in the DES sample.
    for i, SN in enumerate(SNe):
        try:
            # If the source keyword exists
            if (SN['source'] == 'JLA' or SN['source'] == 'SNLS_spec') and SN['name'][2:4] in ['D1', 'D2', 'D3', 'D4']:
                SNe['eval'][i] = True
            elif (SN['source']== 'SNLS_photo') and (SN['name'][2:4] in ['D1', 'D2', 'D3', 'D4'] or (SN['name'][0:2] in ['D1', 'D2', 'D3', 'D4'])):
                SNe['eval'][i] = True
        except:
            # If the source keyword does not exist
            if SN['name'][0:3]=="DES":
                SNe['eval'][i] = True

    print list(SNe['eval']).count(True)
    # Work out which redshift bin each SNe belongs to
    # In numpy.digitize, the bin number starts at 1, so we subtract 1 -- need to check...
    SNe['bin'] = numpy.digitize(SNe['zhel'], z_bin)-1 

    # Build the covariance matrix
    C_nonIa = numpy.zeros(nSNe*3*nSNe*3).reshape(nSNe*3, nSNe*3)

    # It only computes the covariance for the spectroscopically confirmed SNLS SNe
    # We assume that covariance between redshift bins is uncorrelated
    # Within a redshift bin, we assume 100% covariance between SNe in that bin

    for i in range(nSNe):
        bin1 = SNe['bin'][i]
        if SNe['eval'][i]:
            print SNe['zhel'][i], bin1, raw_bias[bin1], f_star[bin1], i
        for j in range(nSNe):
            bin2 = SNe['bin'][j]
            if SNe['eval'][j] and SNe['eval'][i] and bin1 == bin2:
                C_nonIa[3*i, 3*j] = (raw_bias[bin1] * f_star[bin1])**2

    # print SNe['bin'][:239]
    # I am unable to reproduce this JLA covariance matrix

    date = JLA.get_date()

    fits.writeto('C_nonIa_%s.fits' % date, numpy.array(C_nonIa), clobber=True)

    return

if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA parameters")

    PARSER.add_option("-s", "--SNlist", dest="SNlist",
                      help="List of SN")

    PARSER.add_option("-l", "--lcfits", dest="lcfits", default="lightCurveFits",
                      help="Key in config file pointing to lightcurve fit parameters")
    
    PARSER.add_option("-j", "--jla", dest="jla", default=False,
                      action='store_true',
                      help="Only use the SNe from the JLA sample")

    (OPTIONS, ARGS) = PARSER.parse_args()

    compute_nonIa(OPTIONS)
