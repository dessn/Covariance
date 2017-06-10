# pylint: disable=E1121,E1101

"""Python program to compute C_host
"""

from optparse import OptionParser
import numpy as np
import astropy.io.fits as fits
import JLA_library as JLA
from astropy.table import Table

class HostCorrection(object):
    """ functions to contruct covariance matrix for uncertainty in host mass correction """

    def __init__(self):
        self.Delta_M = 0.08

    def H_low(self, SNe):
        """ indicator function for SNe with borderline low mass """
        H = np.zeros((3*len(SNe), 1))
        selection = np.where((SNe['3rdvar'] >= 9) & (SNe['3rdvar'] < 10))[0]
        H[3*selection] = 1.
        return H

    def H_high(self, SNe):
        """ indicator function for SNe with borderline high mass """
        H = np.zeros((3*len(SNe), 1))
        selection = np.where((SNe['3rdvar'] > 10) & (SNe['3rdvar'] <= 11))[0]
        H[3*selection] = 1.
        return H

    def covmat_host(self, SNe):
        """ build covariance matrix, including
        diagonal term for SNe with host galaxy mass uncertainty over the boundary
        and one term each for H_high and H_low
        """
        nSNe = len(SNe)
        chost = np.zeros((3*nSNe, 3*nSNe))
        for i in range(nSNe):
            if np.sign(SNe['3rdvar'][i]-SNe['d3rdvar'][i]-10) != np.sign(SNe['3rdvar'][i]
                                                                        +SNe['d3rdvar'][i]-10):
                chost[3*i][3*i] += self.Delta_M**2
        chost += self.Delta_M**2 * np.dot(self.H_low(SNe), self.H_low(SNe).T)
        chost += self.Delta_M**2 * np.dot(self.H_high(SNe), self.H_high(SNe).T)

        return chost

def compute_Chost(options):

    # -----------  Read in the configuration file ------------
    params = JLA.build_dictionary(options.config)
    
    # ----------  Read in the SN light curve fits ------------
    # We get the host mass from the light curve fits, listed in the column called '3rdvar'.
    lcfile = JLA.get_full_path(params[options.lcfits])
    SN_data = Table.read(lcfile, format='fits')

    # ---------- Read in the SNe list -------------------------
    SN_list_long = np.genfromtxt(options.SNlist, usecols=(0), dtype='S30')
    SN_list = [name.replace('lc-', '').replace('_smp', '').replace('.list', '') for name in SN_list_long]
    SN_indices = JLA.reindex_SNe(SN_list, SN_data)
    SN_data = SN_data[SN_indices]

    host_correction = HostCorrection()

    C_host = host_correction.covmat_host(SN_data)

    date = JLA.get_date()

    fits.writeto('C_host_%s.fits' % date, np.array(C_host), clobber=True)

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containing the location of various JLA parameters")

    parser.add_option("-s", "--SNlist", dest="SNlist",
                      help="List of SNe")

    parser.add_option("-l", "--lcfits", dest="lcfits", default="lightCurveFits",
                      help="Key in config file pointing to lightcurve fit parameters")

    (options, args) = parser.parse_args()

    compute_Chost(options)

