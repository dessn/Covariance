"""
Python program to compute the contribution that each uncertainty makes to the
global uncertainty
Currently, we do this for the JLA sample only
See section 6.2 of B14 for a description of the technique
"""

from optparse import OptionParser

def add_covar_matrices(covmatrices,diag):
    """
    Python program that adds the individual covariance matrices into a single matrix
    """

    import numpy
    import astropy.io.fits as fits
    import JLA_library as JLA

    # Read in the terms that account for uncertainties in perculiar velocities, 
    # instrinsic dispersion and, lensing

    # Read in the covariance matrices
    matrices = []
    for matrix in covmatrices:
        matrices.append(fits.getdata(JLA.get_full_path(covmatrices[matrix]), 0))
        # Test for NaNs and replace them with zero
        if numpy.isnan(matrices[-1]).any():
            print 'Found a NaN in %s ... replacing them with zero' % (covmatrices[matrix])
            print numpy.isnan(matrices[-1]).sum()
            matrices[-1][numpy.isnan(matrices[-1])]=0.0

    # Add the matrices
    size = matrices[0].shape
    add = numpy.zeros(size[0]**2.).reshape(size[0], size[0])
    for matrix in matrices:
        add += matrix

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

    sigma = numpy.genfromtxt(JLA.get_full_path(diag),
                             comments='#',
                             usecols=(0, 1, 2),
                             dtype='f8,f8,f8',
                             names=['sigma_coh', 'sigma_lens', 'sigma_pecvel'])

    for i in range(nSNe):
        cov[i, i] += sigma['sigma_coh'][i]**2 + \
        sigma['sigma_lens'][i]**2 + \
        sigma['sigma_pecvel'][i]**2

    return cov

def compute_rel_size(options):
    import numpy
    import astropy.io.fits as fits
    from astropy.table import Table
    import JLA_library as JLA
    from astropy.cosmology import FlatwCDM
    import os
    
    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)

    # ---------- Read in the SNe list -------------------------

    SNeList=numpy.genfromtxt(options.SNlist,usecols=(0,2),dtype='S30,S200',names=['id','lc'])

    for i,SN in enumerate(SNeList):
        SNeList['id'][i]=SNeList['id'][i].replace('lc-','').replace('.list','')

    # -----------  Read in the data JLA --------------------------

    lcfile = JLA.get_full_path(params[options.lcfits])
    SNe = Table.read(lcfile, format='fits')

    nSNe=len(SNe)
    print 'There are %d SNe in this sample' % (nSNe)

    # sort it to match the listing in options.SNlist
    indices = JLA.reindex_SNe(SNeList['id'], SNe)        
    SNe=SNe[indices]

    # ---------- Compute the Jacobian ----------------------
    # The Jacobian is an m by 4 matrix, where m is the number of SNe
    # The columns are ordered in terms of Om, w, alpha and beta

    J=[]
    JLA_result={'Om':0.303,'w':-1.00,'alpha':0.141,'beta':3.102,'M_B':-19.05}
    offset={'Om':0.01,'w':0.01,'alpha':0.01,'beta':0.01,'M_B':0.01}
    nFit=4

    cosmo1 = FlatwCDM(name='SNLS3+WMAP7', H0=70.0, Om0=JLA_result['Om'], w0=JLA_result['w'])

    # Varying Om
    cosmo2 = FlatwCDM(name='SNLS3+WMAP7', H0=70.0, Om0=JLA_result['Om']+offset['Om'], w0=JLA_result['w'])
    J.append(5*numpy.log10((cosmo1.luminosity_distance(SNe['zcmb'])/cosmo2.luminosity_distance(SNe['zcmb']))[:,0]))

    # varying alpha
    J.append(1.0*offset['alpha']*SNe['x1'][:,0])

    # varying beta
    J.append(-1.0*offset['beta']*SNe['color'][:,0])

    # varying M_B

    J.append(offset['M_B']*numpy.ones(nSNe))
    
    J = numpy.matrix(numpy.concatenate((J)).reshape(nSNe,nFit,order='F') * 100.)

    # Set up the covariance matrices

    systematic_terms = ['bias', 'cal', 'host', 'dust', 'model', 'nonia', 'pecvel', 'stat']

    covmatrices = {'bias':params['bias'],
                   'cal':params['cal'],
                   'host':params['host'],
                   'dust':params['dust'],
                   'model':params['model'],
                   'nonia':params['nonia'],
                   'pecvel':params['pecvel'],
                   'stat':params['stat']}


    if options.type in systematic_terms:
        print "Using %s for the %s term" % (options.name,options.type) 
        covmatrices[options.type]=options.name

    # Combine the matrices to compute the full covariance matrix, and compute its inverse
    if options.all:
        #read in the user provided matrix, otherwise compute it, and write it out
        C=fits.getdata(JLA.get_full_path(params['all']))
    else:
        C=add_covar_matrices(covmatrices,params['diag'])
        date=JLA.get_date()
        fits.writeto('C_total_%s.fits' % (date), C, clobber=True)

    Cinv=numpy.matrix(C).I


    # Construct eta, a 3n vector

    eta=numpy.zeros(3*nSNe)
    for i,SN in enumerate(SNe):
        eta[3*i]=SN['mb']
        eta[3*i+1]=SN['x1']
        eta[3*i+2]=SN['color']

    # Construct A, a n x 3n matrix
    A=numpy.zeros(nSNe*3*nSNe).reshape(nSNe,3*nSNe)

    for i in range(nSNe):
        A[i,3*i]=1.0
        A[i,3*i+1]=JLA_result['alpha']
        A[i,3*i+2]=-JLA_result['beta']

    # ---------- Compute W  ----------------------
    # W has shape m * 3n, where m is the number of fit paramaters.

    W=(J.T * Cinv * J).I * J.T* Cinv* numpy.matrix(A)

    # Note that (J.T * Cinv * J) is a m x m matrix, where m is the number of fit parameters

    # ----------- Compute V_x, where x represents the systematic uncertainty

    result=[]

    for term in systematic_terms:
        cov=numpy.matrix(fits.getdata(JLA.get_full_path(covmatrices[term])))
        if 'C_stat' in covmatrices[term]:
            # Add diagonal term from Eq. 13 to the magnitude
            sigma = numpy.genfromtxt(JLA.get_full_path(params['diag']),comments='#',usecols=(0,1,2),dtype='f8,f8,f8',names=['sigma_coh','sigma_lens','sigma_pecvel'])
            for i in range(nSNe):
                cov[3*i,3*i] += sigma['sigma_coh'][i] ** 2 + sigma['sigma_lens'][i] ** 2 + sigma['sigma_pecvel'][i] ** 2



        V=W * cov * W.T
        result.append(V[0,0])

    print '%20s\t%5s\t%5s\t%s' % ('Term','sigma','var','Percentage')
    for i,term in enumerate(systematic_terms):
        if options.type!=None and term==options.type:
            print '* %18s\t%5.4f\t%5.4f\t%4.1f' % (term,numpy.sqrt(result[i]),result[i],result[i]/numpy.sum(result)*100.)
        else:
            print '%20s\t%5.4f\t%5.4f\t%4.1f' % (term,numpy.sqrt(result[i]),result[i],result[i]/numpy.sum(result)*100.)

    print '%20s\t%5.4f' % ('Total',numpy.sqrt(numpy.sum(result)))

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA parameters")

    parser.add_option("-s", "--SNlist", dest="SNlist",
                      help="List of SN")
    
    parser.add_option("-l", "--lcfits", dest="lcfits", default="lightCurveFits",
                      help="Key in config file pointing to lightcurve fit parameters")

    parser.add_option("-a", "--all", dest="all", default=False,
                      action='store_true',
                      help="Provide the summed matrix")

    parser.add_option("-t", "--type", dest="type", default=None,
                      help="Type of systematic to change")

    parser.add_option("-n", "--name", dest="name", default=None,
                      help="Name of the new matrix")



    (options, args) = parser.parse_args()

    compute_rel_size(options)

