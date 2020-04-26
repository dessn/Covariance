"""A python script to compare SALT surfaces
"""
from optparse import OptionParser
import os
import JLA_library as JLA
import numpy as numpy
from astropy.table import Table
from scipy import interpolate

# See wavelength.h
wave_B=4302.57
wave_V=5428.55

def reduced_lambda(x):
    return (x-wave_B) / (wave_V - wave_B)

def CCM(wave,R_V=3.1):
    val=[]
    for w in wave:
        x=10000.0/w
        if x < 1.1 and x > 0.3:
            a=0.574*x**1.61
            b=-0.527*x**1.61
        elif x< 3.3 and x > 1.1:
            y=x-1.82
            a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
            b =     1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
        else:
            # We drop the other terms as we do not reach that far in the UV
            a = 1.752 - 0.316*x - 0.104 / ((x-4.67)**2 + 0.341)
            b =-3.091 + 1.825*x + 1.206 / ((x-4.62)**2 + 0.263)

        val.append(a + b / R_V)
    return numpy.array(val)

#def Fitz99(wave):
    # For R_V=3.1
    # See http://iopscience.iop.org/article/10.1086/316293/pdf
#    from scipy.interpolate import CubicSpline
#    knots=numpy.array([26500, 12200, 6000, 5470, 4670, 4110, 2700, 2600])
#    values=numpy.array([0.265, 0.829, 2.688, 3.055, 3.806, 4.315, 6.265, 6.591])
#    cs = CubicSpline(knots[::-1], values[::-1])       
#    return cs(wave)

def Fitz99_Spline(wave,R_V):
    # See http://iopscience.iop.org/article/10.1086/316293/pdf
    # for lambda > 2700 Angstroms
    from scipy.interpolate import CubicSpline
    knots=numpy.array([26500., 12200., 6000., 5470., 4670., 4110., 2700., 2600.])
    # Ojo!
    # There is a sign error in the last row of Table 4
    values=numpy.array([0.265, 0.829, -0.426+1.0044*R_V, -0.050+1.0016*R_V, 0.701+1.0016*R_V, 1.208+1.0032*R_V-0.00033*R_V**2, 6.265, 6.591])
    # R_V=3.1
    #values=numpy.array([0.265, 0.829, 2.688, 3.055, 3.806, 4.315, 6.265, 6.591])
    cs = CubicSpline(knots[::-1], values[::-1])
    return cs(wave)


def Fitz99(wave,R_V=3.1):
    # See http://iopscience.iop.org/article/10.1086/316293/pdf
    # and earlier references
    x0=4.596
    gam=0.99
    c2=-0.824+4.717/R_V
    c1=2.030-3.007*c2
    c3=3.23
    c4=0.41
    val=[]
          
    for w in wave:
        x=10000./w
        D=x**2 / ((x**2-x0**2)**2+gam**2*x**2)
        if x > 5.9:
            F=0.5392*(x-5.9)**2 + 0.0564*(x-5.9)**3
        else:
            F=0.0
        val.append(c1+c2*x+c3*D+c4*F)
    
    
    return numpy.array(val)


class readSALTsurface:
    """A SALT surface"""
    def __init__(self,f, salt3Model):
        if salt3Model:
            prefix='salt3'
        else:
            prefix='salt2'

        template_0=Table(numpy.genfromtxt(f+prefix+'_template_0.dat',dtype=[('phase',float),
                                                                           ('wave',float),
                                                                           ('flux',float)]))

        template_1=Table(numpy.genfromtxt(f+prefix+'_template_1.dat',dtype=[('phase',float),
                                                                           ('wave',float),
                                                                           ('flux',float)]))

        salt2_lc_relative_variance_0=Table(numpy.genfromtxt(f+prefix+'_lc_relative_variance_0.dat',dtype=[('phase',float),
                                                                           ('wave',float),
                                                                           ('variance',float)]))

        salt2_lc_relative_variance_1=Table(numpy.genfromtxt(f+prefix+'_lc_relative_variance_1.dat',dtype=[('phase',float),
                                                                           ('wave',float),
                                                                           ('variance',float)]))

        salt2_lc_relative_covariance_01=Table(numpy.genfromtxt(f+prefix+'_lc_relative_covariance_01.dat',dtype=[('phase',float),
                                                                           ('wave',float),
                                                                           ('covariance',float)]))


        self.colour_law=numpy.genfromtxt(f+prefix+'_color_correction.dat',dtype=[('coeff',float)],
                                    skip_header=1,skip_footer=3)
        

        self.colour_law_error=numpy.genfromtxt(f+prefix+'_color_dispersion.dat',dtype=[('wave',float),
                                                                                ('sig',float)], 
                                         skip_header=3)
        

        self.surfaceName=f
        self.template_0={}
        self.template_1={}
        self.salt2_lc_relative_variance_0={}
        self.salt2_lc_relative_variance_1={}
        self.salt2_lc_relative_covariance_01={}
        self.colourlaw=[]

        for phase in numpy.unique(template_0['phase']):
            selection=(template_0['phase']==phase)
            self.template_0['%s' % phase]=template_0[selection]['wave','flux']

        for phase in numpy.unique(template_1['phase']):
            selection=(template_1['phase']==phase)
            self.template_1['%s' % phase]=template_1[selection]['wave','flux']

        for phase in numpy.unique(salt2_lc_relative_variance_0['phase']):
            selection=(salt2_lc_relative_variance_0['phase']==phase)
            self.salt2_lc_relative_variance_0['%s' % phase]=salt2_lc_relative_variance_0[selection]['wave','variance']

        for phase in numpy.unique(salt2_lc_relative_variance_1['phase']):
            selection=(salt2_lc_relative_variance_1['phase']==phase)
            self.salt2_lc_relative_variance_1['%s' % phase]=salt2_lc_relative_variance_1[selection]['wave','variance']

        for phase in numpy.unique(salt2_lc_relative_covariance_01['phase']):
            selection=(salt2_lc_relative_covariance_01['phase']==phase)
            self.salt2_lc_relative_covariance_01['%s' % phase]=salt2_lc_relative_covariance_01[selection]['wave','covariance']

        return


def derivative(alpha,surface,reduced_wave):
    d=alpha
    for exponent,coeff in enumerate(surface.colour_law['coeff']):
        d+=(exponent+2)*coeff*reduced_wave**(exponent+1)
    return d
    
def colourLaw(alpha,surface,reduced_wave):
    d=alpha*reduced_wave
    for exponent,coeff in enumerate(surface.colour_law['coeff']):
        d+=coeff*reduced_wave**(exponent+2)
    return d


def compareSALTsurfaces(surface):

    import matplotlib.pyplot as plt

    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)
    
    # -----------  Read in the SALT models -------------------

    # If the SALT model is SALT3 - only for the second surface
    
    surface1=readSALTsurface(JLA.get_full_path(params['model1']),False)
    surface2=readSALTsurface(JLA.get_full_path(params['model2']),options.salt3)

    name="%s-%s" % (surface1.surfaceName.split('/')[-3],surface2.surfaceName.split('/')[-3])

    #print surface1.surfaceName.split('/')[-3]
    #print surface2.surfaceName.split('/')[-3]

    #print surface1.colour_law
    #print surface2.colour_law

    
    # -----------  Plot the surfaces ----------------------
    fig1=plt.figure()
    for axes,x1 in enumerate([-2,0,2]):
        ax1=fig1.add_subplot(3,1,axes+1)
        flux1=surface1.template_0[options.phase]['flux'] + x1 * surface1.template_1[options.phase]['flux']
        flux2=surface2.template_0[options.phase]['flux'] + x1 * surface2.template_1[options.phase]['flux']
        ax1.plot(surface1.template_0[options.phase]['wave'],flux1,color='b')
        ax1.plot(surface2.template_0[options.phase]['wave'],flux2,color='k')
        ax1.text(7000,0.3,"C=0 x1=%2d" % x1)
        ax1b = ax1.twinx()
        # Accounting for the fact that some models are sampled at twice the rate
        if len(surface1.template_0[options.phase]['wave']) < len(surface2.template_0[options.phase]['wave']):
            #flux1_int=interpolate.interp1d(surface1.template_0[options.phase]['wave'],flux1)(surface2.template_0[options.phase]['wave'])
            #ax1b.plot(surface2.template_0[options.phase]['wave'],(flux1_int-flux2)/flux1_int*100.,color='r')
            #ax1.plot(surface2.template_0[options.phase]['wave'],(flux1_int-flux2),color='g')
            ax1b.plot(surface1.template_0[options.phase]['wave'],(flux1-flux2[0::2])/flux1*100.,color='r')
            ax1.plot(surface1.template_0[options.phase]['wave'],(flux1-flux2[0::2]),color='g')

        else:
            ax1b.plot(surface1.template_0[options.phase]['wave'],(flux1-flux2)/flux1*100.,color='r')
            ax1.plot(surface1.template_0[options.phase]['wave'],(flux1-flux2),color='g')

        ax1b.plot([2000,9200],[0,0],'m--',marker=None)
        ax1b.set_ylim(-20,20)
        ax1b.tick_params(axis='y', labelcolor='r')
        
    ax1.set_xlabel("wavelength ($\AA$)")
    ax1.set_ylabel("flux density")
    ax1b.set_ylabel('% Difference', color='r')
        
    plt.savefig("%s_SED.png" % (name))

       
    # -----------  Plot the colour laws ----------------------

    # See salt2extinction.cc
    # Note the extrapolation

#    /*
#    ========================================================
#    VERSION 1
#    ========================================================
#    if(l_B<=l<=l_R)
#    ext = exp( color * constant * ( alpha*l + params(0)*l^2 + params(1)*l^3 + ... ))
#    = exp( color * constant * P(l) )
#    alpha = 1-params(0)-params(1)-...
#    if(l>l_R)
#    ext = exp( color * constant * ( P(l_R) + P'(l_R)*(l-l_R) ) )
#    if(l<l_B)
#    ext = exp( color * constant * ( P(l_B) + P'(l_B)*(l-l_B) ) )
#  
#    ======================================================== 
#    */

    constant=0.4 * numpy.log(10)
    fig3=plt.figure()
    ax3=fig3.add_subplot(111)
    wave1=surface1.template_0[options.phase]['wave']
    wave2=surface2.template_0[options.phase]['wave']
    wave1_min=2800.
    wave1_max=7000.
    wave2_min=2800.
    wave2_max=7000.
    
    wave1_min_reduced=reduced_lambda(wave1_min)
    wave1_max_reduced=reduced_lambda(wave1_max)
    wave2_min_reduced=reduced_lambda(wave2_min)
    wave2_max_reduced=reduced_lambda(wave2_max)
    
    # See salt2extinction.h
    reduced_wave1=reduced_lambda(wave1)
    reduced_wave2=reduced_lambda(wave2)

    # Model 1
    alpha1=1.0
    
    # There are 4 co-efficients in the colour law
    for coeff in surface1.colour_law['coeff']:
        alpha1-=coeff

    p1=numpy.zeros(len(reduced_wave1))
        
    # Compute derivatives for extrapolations
    p1_derivative_min=derivative(alpha1,surface1,wave1_min_reduced)
    p1_derivative_max=derivative(alpha1,surface1,wave1_max_reduced)

    # Compute colour law at the points of extrapolations
    p1_wave_min_reduced=colourLaw(alpha1,surface1,wave1_min_reduced)
    p1_wave_max_reduced=colourLaw(alpha1,surface1,wave1_max_reduced)

    for index,rl in enumerate(reduced_wave1):
        if rl < wave1_min_reduced:
            p1[index]=p1_wave_min_reduced+p1_derivative_min*(rl-wave1_min_reduced)
        elif rl > wave1_max_reduced:
            p1[index]=p1_wave_max_reduced+p1_derivative_max*(rl-wave1_max_reduced)
        else:
            p1[index]=colourLaw(alpha1,surface1,rl)
    
    # Model 2
    alpha2=1.0

    for coeff in surface2.colour_law['coeff']:
        alpha2-=coeff

    p2=numpy.zeros(len(reduced_wave2))
            
    # Compute derivatives for extrapolations
    p2_derivative_min=derivative(alpha2,surface2,wave2_min_reduced)
    p2_derivative_max=derivative(alpha2,surface2,wave2_max_reduced)

    # Compute colour law at the points of extrapolations
    p2_wave_min_reduced=colourLaw(alpha2,surface2,wave2_min_reduced)
    p2_wave_max_reduced=colourLaw(alpha2,surface2,wave2_max_reduced)

    for index,rl in enumerate(reduced_wave2):
        if rl < wave2_min_reduced:
            p2[index]=p2_wave_min_reduced+p2_derivative_min*(rl-wave2_min_reduced)
        elif rl > wave2_max_reduced:
            p2[index]=p2_wave_max_reduced+p2_derivative_max*(rl-wave2_max_reduced)
        else:
            p2[index]=colourLaw(alpha2,surface2,rl)
    
    # See Fig.3 of B14.
    # p1 and p2 are the log (colour law)
    
    C=-0.1

    # Allowing for the fact that the SEDs might have finer sampling
    if len(p1) > len(surface1.colour_law_error['sig']):
        p1_s2=p1[0::2]
    else:
        p1_s2=p1

    A1_wave=p1_s2*C
    A1_wave_err_plus=(p1_s2+surface1.colour_law_error['sig'])*C
    A1_wave_err_minus=(p1_s2-surface1.colour_law_error['sig'])*C
    ax3.fill_between(surface1.colour_law_error['wave'], A1_wave_err_plus, A1_wave_err_minus, alpha=0.4,label='model1+err')
    ax3.plot(surface1.colour_law_error['wave'], A1_wave,label='model1')
    
    
    A2_wave=p2*C
    ax3.plot(wave2, A2_wave, label='model2')

    # Plot CCM R_V=3.1
    E_BV=0.1
    R_V=3.1
    a_wave=E_BV * R_V * CCM(wave1, R_V)
    a_B=E_BV * R_V * CCM(numpy.array([wave_B]),R_V)
    ax3.plot(wave1,a_wave-a_B,label='CCM R_V=3.1')

    #CCM R_V=1.0
    R_V=1.0
    a_wave=E_BV * R_V * CCM(wave1, R_V)
    a_B=E_BV * R_V * CCM(numpy.array([wave_B]),R_V)
    ax3.plot(wave1,a_wave-a_B,label='CCM R_V=1.0')

    # F99 R_V=3.1
    
    # Fitz99 UV extiction - note the limited wavelength range
    # https://iopscience.iop.org/article/10.1086/316293/pdf
    # http://adsabs.harvard.edu/abs/1990ApJS...72..163FE_BV=0.1
    R_V=3.1
    wave=numpy.arange(2000,2700,1.0)
    ax3.plot(wave, E_BV*(Fitz99(wave,R_V)-1.0),label="Fitz99 UV R_V=%3.1f" % R_V)
    # Fitz99 - Spline function in the optical
    # https://iopscience.iop.org/article/10.1086/316293/pdf
    wave=numpy.arange(2700,9200,1.0)
    a_B=E_BV*(Fitz99_Spline(numpy.array([wave_B]),R_V))
    ax3.plot(wave,E_BV*Fitz99_Spline(wave,R_V)-a_B,label="Fitz99 Spline R_V=%3.1f" % R_V)
    
    ax3.legend()
    ax3.set_xlabel("wavelength ($\AA$)")
    ax3.set_ylim(-0.3,0.8)
    plt.savefig("%s_colourlaw.png" % (name))

    # -----------  Plot examples of the impact of colour ----------------------
    # Assume x1=0
    # Note
    # The colour laws p1 and p2 have the absorption in the B band subtracted
    # The units are magnitudes
    # Are we correctly applyng the colour law?

    fig2=plt.figure()
    for axes,C in enumerate([-0.1,0,0.1]):
        ax2=fig2.add_subplot(3,1,axes+1)
        flux1=surface1.template_0[options.phase]['flux'] * numpy.exp(C*constant*p1)
        flux2=surface2.template_0[options.phase]['flux'] * numpy.exp(C*constant*p2)
        ax2.plot(surface1.template_0[options.phase]['wave'],flux1,color='b')
        ax2.plot(surface2.template_0[options.phase]['wave'],flux2,color='k')
        ax2.text(7000,0.3,"C=%4.1f x1=0.0" % C)
        ax2b = ax2.twinx()
        ax2b.plot([2000,9200],[0,0],'m--',marker=None)
        if len(surface1.template_0[options.phase]['wave']) < len(surface2.template_0[options.phase]['wave']):
#            flux1_int=interpolate.interp1d(surface1.template_0[options.phase]['wave'],flux1)(surface2.template_0[options.phase]['wave'])
#            ax2b.plot(surface2.template_0[options.phase]['wave'],(flux1_int-flux2)/flux1_int*100.,color='r')
#            ax2.plot(surface2.template_0[options.phase]['wave'],(flux1_int-flux2),color='g')
            ax2b.plot(surface1.template_0[options.phase]['wave'],(flux1-flux2[0::2])/flux1*100.,color='r')
            ax2.plot(surface1.template_0[options.phase]['wave'],(flux1-flux2[0::2]),color='g')
        else:
            ax2b.plot(surface1.template_0[options.phase]['wave'],(flux1-flux2)/flux1*100.,color='r')
            ax2.plot(surface1.template_0[options.phase]['wave'],(flux1-flux2),color='g')
        ax2b.set_ylim(-20,20)
        ax2b.tick_params(axis='y', labelcolor='r')

    ax2.set_xlabel("wavelength ($\AA$)")
    ax2.set_ylabel("flux density")
    ax2b.set_ylabel('% Difference', color='r')

    plt.savefig("%s_colour_SED.png" % (name))

    # Examine the differenece in the variance 
    fig4=plt.figure()
    ax4=fig4.add_subplot(3,1,1)
    ax4.plot(surface1.salt2_lc_relative_variance_0[options.phase2]['wave'],surface1.salt2_lc_relative_variance_0[options.phase2]['variance'],color='b')
    ax4.plot(surface2.salt2_lc_relative_variance_0[options.phase2]['wave'],surface2.salt2_lc_relative_variance_0[options.phase2]['variance'],color='k')
    ax4b = ax4.twinx()
    diff=surface1.salt2_lc_relative_variance_0[options.phase2]['variance']-surface2.salt2_lc_relative_variance_0[options.phase2]['variance']
    ax4b.plot(surface1.salt2_lc_relative_variance_0[options.phase2]['wave'],diff/surface1.salt2_lc_relative_variance_0[options.phase2]['variance']*100,color='r')
    ax4b.set_ylim(-20,20)
    ax4b.tick_params(axis='y', labelcolor='r')

    ax4.set_xlabel("wavelength ($\AA$)")
    ax4.set_ylabel("var (comp 0)")
    ax4b.set_ylabel('% Difference', color='r')

    ax4=fig4.add_subplot(3,1,2)
    ax4.plot(surface1.salt2_lc_relative_variance_1[options.phase2]['wave'],surface1.salt2_lc_relative_variance_1[options.phase2]['variance'],color='b')
    ax4.plot(surface2.salt2_lc_relative_variance_1[options.phase2]['wave'],surface2.salt2_lc_relative_variance_1[options.phase2]['variance'],color='k')
    ax4b = ax4.twinx()
    diff=surface1.salt2_lc_relative_variance_1[options.phase2]['variance']-surface2.salt2_lc_relative_variance_1[options.phase2]['variance']
    ax4b.plot(surface1.salt2_lc_relative_variance_1[options.phase2]['wave'],diff/surface1.salt2_lc_relative_variance_1[options.phase2]['variance']*100,color='r')
    ax4b.set_ylim(-20,20)
    ax4b.tick_params(axis='y', labelcolor='r')

    ax4.set_xlabel("wavelength ($\AA$)")
    ax4.set_ylabel("var (comp 1)")
    ax4b.set_ylabel('% Difference', color='r')

    ax4=fig4.add_subplot(3,1,3)
    ax4.plot(surface1.salt2_lc_relative_covariance_01[options.phase2]['wave'],surface1.salt2_lc_relative_covariance_01[options.phase2]['covariance'],color='b')
    ax4.plot(surface2.salt2_lc_relative_covariance_01[options.phase2]['wave'],surface2.salt2_lc_relative_covariance_01[options.phase2]['covariance'],color='k')
    ax4b = ax4.twinx()
    diff=surface1.salt2_lc_relative_covariance_01[options.phase2]['covariance']-surface2.salt2_lc_relative_covariance_01[options.phase2]['covariance']
    ax4b.plot(surface1.salt2_lc_relative_covariance_01[options.phase2]['wave'],diff/surface1.salt2_lc_relative_covariance_01[options.phase2]['covariance']*100,color='r')
    ax4b.set_ylim(-20,20)
    ax4b.tick_params(axis='y', labelcolor='r')

    ax4.set_xlabel("wavelength ($\AA$)")
    ax4.set_ylabel("covar (comp 01)")
    ax4b.set_ylabel('% Difference', color='r')

    plt.savefig("%s_lc_variance.png" % (name))


    plt.show()
    plt.close()

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="SALTmodels.config",
                      help="Directories containing the SALT models")

    parser.add_option("-p", "--phase", dest="phase", default=0.0,
                      help="Lightcurve phase")

    parser.add_option("--phase2", dest="phase2", default=0.5,
                      help="Lightcurve phase")

    parser.add_option("-3", "--salt3", dest="salt3", default=False, action="store_true",
                      help="SALT3 model")

    (options, args) = parser.parse_args()


    compareSALTsurfaces(options)
