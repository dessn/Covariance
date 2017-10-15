# pylint: disable=E1121,E1101

"""Python program to compute C_pecvel
need to specify the list of SNe and
lightcurve fit parameters (if not subset of JLA)
"""

from optparse import OptionParser
import numpy as np
import astropy.io.fits as fits
import JLA_library as JLA

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

def sin_d(x):
    """ sine in degrees for convenience """
    return np.sin(x*np.pi/180.)

def cos_d(x):
    """ cosine in degrees, for convenience """
    return np.cos(x*np.pi/180.)

def dmdz(z):
    """ converts uncertainty in redshift to uncertainty in magnitude, for clarity"""
    return 5./(np.log(10)*z)

class VelocityCorrection(object):
    """ functions to correct for peculiar velocities using 2M++ velocity field (Carrick+ 2015) """

    def __init__(self, velfield):
        """ sets up constants and velocity field; includes:
        params of CMB dipole w.r.t. heliocentric frame (Bennett+ 2003)
        params of LG w.r.t. heliocentric frame (Tully+ 2008 table 3)
        all velocities in km/s; angles in degrees
        sets up 2M++ velocity field is in galactic cartesian coordinates,
        in local group frame and in Mpc/h
        this is linear with coefficient beta
        we perturb beta to measure the effect to compute the covariance matrix
        """

        self.c = 3.00*10**5

        self.v_helio = 371.
        self.l_h = 263.85 # galactic longitude
        self.b_h = 48.25 # galactic latitude

        self.v_LG = 318.
        self.l_LG = 105.7467
        self.b_LG = -5.9451

        self.velocity_field = np.load(velfield)
        self.beta = 0.43
        self.beta_err = 0.1
        self.r_plus = 1 + self.beta_err/self.beta
        self.r_minus = 1 - self.beta_err/self.beta

    def convert_helio_to_LG(self, v, l, b):
        """ converts velocity from heliocentric frame to Local Group frame """

        return v - self.v_LG*(sin_d(b)*sin_d(self.b_LG)
                              + cos_d(b)*cos_d(self.b_LG)*cos_d(l-self.l_LG))

    def lookup_velocity(self, z_h, l, b):
        """ look up velocity field to find peculiar velocity vector
        input needs to be heliocentric position (and in degrees)
        output is vector in galactic cartesian coordinates, in CMB frame
        """

        cz_LG = self.convert_helio_to_LG(self.c*z_h, l, b)
        x = cz_LG * cos_d(l) * cos_d(b)
        y = cz_LG * sin_d(l) * cos_d(b)
        z = cz_LG * sin_d(b)
        x, y, z = x/100, y/100, z/100 # to convert from km/s to Mpc/h
        i = int(round(128. + 256./400*x))
        j = int(round(128. + 256./400*y))
        k = int(round(128. + 256./400*z))
        try:
            vpec = self.velocity_field[:, i, j, k]
        except IndexError:
            vpec = np.array([0, 0, 0])#outside velocity field; approximate as zero
        return vpec

    def correct_redshift(self, z_h, vpec, l, b):
        """ convert helio redshift to cosmological redshift (zbar; in CMB frame)
        input needs to be vector in galactic cartesian coordinates, w.r.t. CMB frame
        components are: heliocentric motion, peculiar velocity (in radial direction)
        """
        helio_corr = self.v_helio/self.c*((sin_d(b)*sin_d(self.b_h)
                                           + cos_d(b)*cos_d(self.b_h)*cos_d(l-self.l_h)))
        pec_corr = vpec.dot(np.array([cos_d(l)*cos_d(b),
                                      sin_d(l)*cos_d(b), sin_d(b)]))/self.c
        corr_term = 1 - helio_corr + pec_corr
        return (1+z_h)/corr_term - 1

    def apply(self, SNe, to_write=False):
        """ correction to SNe"""
        z_value, z_err = [], []

        for sn in SNe:
            c = SkyCoord(sn['ra']*u.degree, sn['dec']*u.degree)
            gc = c.galactic
            vpec = self.lookup_velocity(sn['zhel'], gc.l.degree, gc.b.degree)
            z_c = self.correct_redshift(sn['zhel'], vpec, gc.l.degree, gc.b.degree)
            z_plus = self.correct_redshift(sn['zhel'], self.r_plus*vpec,
                                           gc.l.degree, gc.b.degree)
            z_minus = self.correct_redshift(sn['zhel'], self.r_minus*vpec,
                                            gc.l.degree, gc.b.degree)
            if JLA.survey(sn) == 'nearby':
                z_value.append(z_c)
                z_err.append(np.mean([z_plus - z_c, z_c - z_minus]))
            else:
                z_value.append(sn['zcmb'])
                z_err.append(0)
                
        if to_write:
            np.savetxt('z_CMB_corrected.txt', np.array(z_value))
        return z_value, z_err


    def covmat_pecvel(self, SNe):
        """ build covariance matrix """
        z, z_err = self.apply(SNe)
        nSNe = len(SNe)
        cpecvel = np.zeros((3*nSNe, 3*nSNe))

        for i in range(nSNe):
            for j in range(i, nSNe):
                cpecvel[3*i, 3*j] = dmdz(z[i])*dmdz(z[j]) * z_err[i] * z_err[j]
                cpecvel[3*j, 3*i] = cpecvel[3*i, 3*j]

        return cpecvel

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containing the location of various JLA parameters")

    parser.add_option("-s", "--SNlist", dest="SNlist",
                      help="List of SNe")

    parser.add_option("-l", "--lcfits", dest="lcfits", default="lightCurveFits",
                      help="Key in config file pointing to lightcurve fit parameters")

    (options, args) = parser.parse_args()

    params = JLA.build_dictionary(options.config)

    lcfile = JLA.get_full_path(params[options.lcfits])
    SN_data = Table.read(lcfile, format='fits')

    SN_list_long = np.genfromtxt(options.SNlist, usecols=(0), dtype='S30')
    SN_list = [name.replace('lc-', '').replace('.list', '').replace('_smp','') for name in SN_list_long] 
  
    SN_indices = JLA.reindex_SNe(SN_list, SN_data)
    SN_data = SN_data[SN_indices]

    velfile = JLA.get_full_path(params['velocityField'])
    vel_correction = VelocityCorrection(velfile)
    #z_correction = vel_correction.apply(SN_data)

    C_pecvel = vel_correction.covmat_pecvel(SN_data)

    date = JLA.get_date()

    fits.writeto('C_pecvel_%s.fits' % date, np.array(C_pecvel), clobber=True)
