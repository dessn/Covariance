# pylint: disable=E1121,E1101

"""Python program to compute C_pecvel for the low-z SNe in JLA sample
"""

import os
import time
from optparse import OptionParser
import numpy as np
import pyfits
import JLA_library as JLA

from astropysics.coords import ICRSCoordinates, GalacticCoordinates

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
    """ functions to do stuff ?? with 2M++ velocity field (Carrick+ 2015) """

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

        cz_LG = self.convert_helio_to_LG(self, self.c*z_h, l, b)
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
            print 'Warning - outside velocity field, retreating to edge ', k
            k = min(256, k) # approximate at edge
            vpec = self.velocity_field[:, i, j, k]

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
            coords = ICRSCoordinates(sn['ra'], sn['dec'])
            gcoords = coords.convert(GalacticCoordinates)
            vpec = self.lookup_velocity(self, sn['z_helio'], gcoords.l.d, gcoords.b.d)
            z_c = self.correct_redshift(self, sn['z_helio'], vpec, gcoords.l.d, gcoords.b.d)
            z_value.append(z_c)
            z_plus = self.correct_redshift(self, sn['z_helio'], self.r_plus*vpec,
                                           gcoords.l.d, gcoords.b.d)
            z_minus = self.correct_redshift(self, sn['z_helio'], self.r_minus*vpec,
                                            gcoords.l.d, gcoords.b.d)
            z_err.append(np.mean([z_plus - z_c, z_c - z_minus]))

        if to_write:
            np.savetxt('z_CMB_corrected_%s.txt', np.array(z_value))
        return z_value, z_err


    def covmat_pecvel(self, SNe):
        """ build covariance matrix """
        z, z_err = self.apply(SNe)
        nSNe = len(SNe)
        cpecvel = np.zeros((3*nSNe, 3*nSNe))

        for i in range(nSNe):
            for j in range(i, nSNe):
                cpecvel[3*i, 3*j] = dmdz(z[i-616])*dmdz(z[j-616]) * z_err[i-616] * z_err[j-616]
                cpecvel[3*j, 3*i] = cpecvel[3*i, 3*j]

        return cpecvel

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containing the location of various JLA parameters")

    parser.add_option("-l", "--lcdir", dest="lcdir", default="use_JLA",
                      help="File containing SN light curve fits (formatted as in jla_lcparams.txt)")

    (options, args) = parser.parse_args()

    params = JLA.build_dictionary(options.config)
    JLAHOME = os.environ["JLA"]

    if options.lcdir == "use_JLA":
        lcfits = JLAHOME+'/'+params['lightCurveFits']
    else:
        lcfits = options.lcdir

    SN_data = np.genfromtxt(lcfits, skip_header=1, usecols=(0, 1, 2, 18, 19),
                            dtype='S10, f8, f8, f8, f8',
                            names=['id', 'z_CMB', 'z_helio', 'ra', 'dec'])

    if options.lcdir == "use_JLA":
        SN_data = SN_data[616:] # to use only low-z SNe in JLA only

    velfld = JLAHOME+'/'+params['velocityField']
    vel_correction = VelocityCorrection(velfld)
    #z_correction = vel_correction.apply(SN_data)

    C_pecvel = vel_correction.covmat_pecvel(SN_data)

    t = time.gmtime(time.time())
    date = '%4d%02d%02d' % (t[0], t[1], t[2])

    pyfits.writeto('C_pecvel_%s.fits' % (date), np.array(C_pecvel), clobber=True)
