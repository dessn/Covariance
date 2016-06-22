"""
Python program that adds the individual covariance matrices into a single matrix
"""

from optparse import OptionParser

def add_covar_matrices(options):
    """
    Python program that adds the individual covariance matrices into a single matrix
    """

    import time
    import numpy
    import astropy.io.fits as fits
    import JLA_library as JLA

    params = JLA.build_dictionary(options.config)

    # Read in the terms that account for uncertainties in perculiar velocities, 
    # instrinsic dispersion and, lensing

    # Read in the covariance matrices
    matrices = []

    systematic_terms = ['bias',  'host', 'dust', 'model', 'nonia', 'pecvel', 'stat']#'cal',

    covmatrices = {'bias':params['bias'],
                   #'cal':params['cal'],
                   'host':params['host'],
                   'dust':params['dust'],
                   'model':params['model'],
                   'nonia':params['nonia'],
                   'pecvel':params['pecvel'],
                   'stat':params['stat']}

    for term in systematic_terms:
        matrices.append(fits.getdata(JLA.get_full_path(covmatrices[term]), 0))

    # Add the matrices
    size = matrices[0].shape
    add = numpy.zeros(size[0]**2.).reshape(size[0], size[0])
    for matrix in matrices:
        add += matrix

    # Write out this matrix. This is C_eta in qe. 13 of B14

    date=JLA.get_date()

    fits.writeto('C_eta_%s.fits' % (date), add, clobber=True)

    # Compute A

    nSNe = size[0]/3

    jla_results = {'Om':0.303, 'w':-1.027, 'alpha':0.141, 'beta':3.102}

    arr = numpy.zeros(nSNe*3*nSNe).reshape(nSNe, 3*nSNe)

    for i in range(nSNe):
        arr[i, 3*i] = 1.0
        arr[i, 3*i+1] = jla_results['alpha']
        arr[i, 3*i+2] = -jla_results['beta']

    cov = numpy.matrix(arr) * numpy.matrix(add) * numpy.matrix(arr).T

    # Add the diagonal terms

    sigma = numpy.genfromtxt(JLA.get_full_path(params['diag']),
                             comments='#',
                             usecols=(0, 1, 2),
                             dtype='f8,f8,f8',
                             names=['sigma_coh', 'sigma_lens', 'sigma_pecvel'])

    for i in range(nSNe):
        cov[i, i] += sigma['sigma_coh'][i]**2 + \
        sigma['sigma_lens'][i]**2 + \
        sigma['sigma_pecvel'][i]**2

    fits.writeto('C_total_%s.fits' % (date), cov, clobber=True)

    return

if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", 
                      default="Shuvo.config",
                      help="Parameter file containting the location of the JLA files")

    (OPTIONS, ARGS) = PARSER.parse_args()

    add_covar_matrices(OPTIONS)

