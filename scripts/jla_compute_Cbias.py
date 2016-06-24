"""
Python program to compute C_bias
"""

from optparse import OptionParser

def poly(x,p):
    val=0
    for i in range(len(p)):
        val+=p[i]*x**i
    return val

def chebyshev(x,p):
    # To generalise, see https://en.wikipedia.org/wiki/Chebyshev_polynomials
    # Only for test purposes
    val=p[0]+p[1]*x+p[2]*(2*x**2-1)
    return val
    

def residuals(p,y,x,e,fn):
    if fn=='poly':
        return (y-poly(x,p))/e
    elif fn=='cheb':
        return (y-chebyshev(x,p))/e


def compute_bias(options):

    import numpy
    import astropy.io.fits as fits
    import JLA_library as JLA
    from astropy.table import Table
    from astropy.cosmology import FlatwCDM
    from  scipy.optimize import leastsq
    import matplotlib.pyplot as plt
    from scipy.stats import t


    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)

    # -----------  Read in the SN ordering ------------------------
    SNeList = Table(numpy.genfromtxt(options.SNlist,
                               usecols=(0, 2),
                               dtype='S30,S200',
                               names=['id', 'lc']))
    nSNe = len(SNeList)

    for i, SN in enumerate(SNeList):
        SNeList['id'][i] = SNeList['id'][i].replace('lc-', '').replace('.list', '')
        
    lightCurveFits = JLA.get_full_path(params['lightCurveFits'])
    SNe = Table.read(lightCurveFits, format='fits')

    print 'There are %d SNe' % (nSNe)

    indices = JLA.reindex_SNe(SNeList['id'], SNe)
    SNe=SNe[indices]
    # Add a column that records the error in the bias
    SNe['e_bias'] = numpy.zeros(nSNe,'f8')


    # Read in the points from B14 figure 
    # Fit a polynomial to the data
    # Determine the uncertainties

    bias = numpy.genfromtxt(JLA.get_full_path(params['biasPolynomial']),
                                  skip_header=3,
                                  usecols=(0, 1, 2, 3),
                                  dtype='S10,f8,f8,f8',
                                  names=['sample', 'redshift', 'bias', 'e_bias'])

    
    if options.plot:
        fig=plt.figure()
        ax=fig.add_subplot(111)
        colour={'nearby':'b','SNLS':'r','SDSS':'g'}

    for sample in numpy.unique(bias['sample']):
        selection=(bias['sample']==sample)
        guess=[0,0,0]
        
        plsq=leastsq(residuals, guess, args=(bias[selection]['bias'],
                                             bias[selection]['redshift'],
                                             bias[selection]['e_bias'],
                                             'poly'), full_output=1)
#
        if plsq[4] in [1,2,3,4]:
            print 'Solution for %s found' % (sample)

        if options.plot:
            ax.errorbar(bias[selection]['redshift'],
                    bias[selection]['bias'],
                    yerr=bias[selection]['e_bias'],
                    ecolor='k',
                    color=colour[sample],
                    fmt='o',
                    label=sample)
            z=numpy.arange(numpy.min(bias[selection]['redshift']),numpy.max(bias[selection]['redshift']),0.001)
            ax.plot(z,poly(z,plsq[0]),color=colour[sample])

        # For each SNe, determine the uncerainty in the correction. We use the covariance martix
        # prediction bounds for the fitted curve. 
        # https://www.astro.rug.nl/software/kapteyn/kmpfittutorial.html
        
        # Compute the chi-sq.
        chisq=(((bias[selection]['bias']-poly(bias[selection]['redshift'],plsq[0]))/bias[selection]['e_bias'])**2.).sum()
        dof=selection.sum()-len(guess)
        print "Reduced chi-square value for sample %s is %5.2e" % (sample, chisq / dof)

        alpha=0.315 # Confidence interval is 100 * (1-alpha)
        # Compute the upper alpha/2 vallue for the student t distribution with dof
        thresh=t.ppf((1-alpha/2.0), dof)
        
        if options.plot:
            # The following is only valid for polynomial fitting functions
            upper_curve=[]
            lower_curve=[]
            for x in z:
                vect=numpy.matrix([1,x,x**2.])
                offset=thresh * numpy.sqrt(chisq / dof * (vect*numpy.matrix(plsq[1])*vect.T)[0,0])
                upper_curve.append(poly(x,plsq[0])+offset)
                lower_curve.append(poly(x,plsq[0])-offset)

            ax.plot(z,lower_curve,'--',color=colour[sample])
            ax.plot(z,upper_curve,'--',color=colour[sample])

        # Compute the error in the bias
        # We increase the absolute vlaue
        # In other words, if the bias is negative, we subtract the error to make it even more negative
        # We assume 100% correlation between SNe
        for i,SN in enumerate(SNe):
            if JLA.survey(SN) == sample:
                redshift = SN['zcmb']
                vect = numpy.matrix([1,redshift,redshift**2.])
                if poly(redshift,plsq[0]) > 0:
                    sign = 1
                else:
                    sign = -1

                SNe['e_bias'][i] = sign * thresh * numpy.sqrt(chisq / dof * (vect*numpy.matrix(plsq[1])*vect.T)[0,0])

    if options.plot:
        ax.legend()
        plt.show()
        plt.close()

    # Compute the bias matrix
    # 

    date = JLA.get_date()
    Zero=numpy.zeros(nSNe)
    H=numpy.concatenate((SNe['e_bias'],Zero,Zero)).reshape(3,nSNe).ravel(order='F')

    C_bias = numpy.matrix(H)

    fits.writeto('C_bias_%s.fits' % (date),C_bias.T*C_bias,clobber=True) 

    return None

if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA parameters")

    PARSER.add_option("-s", "--SNlist", dest="SNlist", 
                      help="List of SN")

    PARSER.add_option("-p", "--plot", dest="plot",default=False,
                      action='store_true',
                      help="List of SN")


    (options, args) = PARSER.parse_args()

    compute_bias(options)
