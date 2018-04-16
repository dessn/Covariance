""" Python program that computes the diagonal terms that come from lensing, peculiar velocities and instrinsic uncertainties
"""

from optparse import OptionParser
import numpy as np
import astropy.io.fits as fits
import JLA_library as JLA
from astropy.table import Table

coh_dict = {'SDSS': 0.108, 'SNLS': 0.08, 'nearby': 0.134, 'high-z': 0.1, 'DES': 0.1}

def compute_diag(SNe):
    lens = 0.055
    #c = 3e5 #km/s
    sigma_lens = 0.055 * SNe['zcmb']
    sigma_pecvel = 5*7.6e-4/(np.log(10) *SNe['zcmb']) # See eq. 13 in B14
    sigma_coh = np.array([coh_dict[JLA.survey(sn)] for sn in SNe])
    return np.column_stack((sigma_coh, sigma_lens, sigma_pecvel))

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containing the location of various JLA parameters")

    parser.add_option("-s", "--SNlist", dest="SNlist",
                      help="List of SNe")

    parser.add_option("-l", "--lcfits", dest="lcfits", default="lightCurveFits",
                      help="Key in config file pointing to lightcurve fit parameters")
    
    parser.add_option("-o", "--output", dest="output",default="sigma_mu.txt", 
                  help="Output")

    (options, args) = parser.parse_args()

    params = JLA.build_dictionary(options.config)
    
    lcfile = JLA.get_full_path(params[options.lcfits])
    SN_data = Table.read(lcfile, format='fits')

    SN_list_long = np.genfromtxt(options.SNlist, usecols=(0), dtype='S30')
    SN_list = [name.replace('lc-', '').replace('.list', '').replace('_smp','') for name in SN_list_long]
    SN_indices = JLA.reindex_SNe(SN_list, SN_data)
    SN_data = SN_data[SN_indices]

    sigma_diag = compute_diag(SN_data)

    np.savetxt(options.output,sigma_diag, header='coh lens pecvel')
