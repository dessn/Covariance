#!/Users/clidman/Science/Programs/Ureka/variants/common/bin/python
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
    import pyfits
    import JLA_library as JLA

    params = JLA.build_dictionary(options.config)

    # Read in the terms that account for uncertainties in perculiar velocities, 
    # instrinsic dispersion and, lensing

    # Read in the covariance matrices
    matrices = []
    for matrix in params["matrices"].split(','):
        matrices.append(pyfits.getdata(JLA.get_full_path(matrix), 0))

    # Add the matrices
    size = matrices[0].shape
    add = numpy.zeros(size[0]**2.).reshape(size[0], size[0])
    for matrix in matrices:
        add += matrix

    # Write out this matrix. This is C_eta in qe. 13 of B14

    ymd = time.gmtime(time.time())
    date = '%4d%02d%02d' % (ymd[0], ymd[1], ymd[2])

    pyfits.writeto('C_eta_%s.fits' % (date), add, clobber=True)

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

    sigma = numpy.genfromtxt(JLA.get_full_path(options.diagonal),
                             comments='#',
                             usecols=(0, 1, 2),
                             dtype='f8,f8,f8',
                             names=['sigma_coh', 'sigma_lens', 'z'])

    sigma_pecvel = (5 * 150 / 3e5) / (numpy.log(10.) * sigma['z'])

    for i in range(nSNe):
        cov[i, i] += sigma['sigma_coh'][i]**2 + \
        sigma['sigma_lens'][i]**2 + \
        sigma_pecvel[i]**2

    pyfits.writeto('C_%s.fits' % (date), cov, clobber=True)

    return

if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", 
                      default="Shuvo.config",
                      help="Parameter file containting the location of the JLA files")

    PARSER.add_option("-d", "--diagonal", dest="diagonal", 
                      default="sigma_mu.txt",
                      help="Diagonal terms")

    (OPTIONS, ARGS) = PARSER.parse_args()

    add_covar_matrices(OPTIONS)

