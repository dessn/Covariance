"""A python script to compare SALT surfaces
"""
from optparse import OptionParser
import os
import JLA_library as JLA
import numpy as numpy
from astropy.table import Table

# See wavelength.h
wave_B=4302.57
wave_V=5428.55

def reduced_lambda(x):
    # What is this ???
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

def Fitz99(wave):
    # For R_V=3.1
    # See http://iopscience.iop.org/article/10.1086/316293/pdf
    from scipy.interpolate import CubicSpline
    knots=numpy.array([26500, 12200, 6000, 5470, 4670, 4110, 2700, 2600])
    values=numpy.array([0.265, 0.829, 2.688, 3.055, 3.806, 4.315, 6.265, 6.591])
    cs = CubicSpline(knots[::-1], values[::-1])       
    return cs(wave)

class readSALTsurface:
    """A SALT surface"""
    def __init__(self,f):
        template_0=Table(numpy.genfromtxt(f+'/salt2_template_0.dat',dtype=[('phase',float),
                                                                           ('wave',float),
                                                                           ('flux',float)]))

        template_1=Table(numpy.genfromtxt(f+'/salt2_template_1.dat',dtype=[('phase',float),
                                                                           ('wave',float),
                                                                           ('flux',float)]))

        self.colour_law=numpy.genfromtxt(f+'/salt2_color_correction.dat',dtype=[('coeff',float)],
                                    skip_header=1,skip_footer=3)
        
        self.template_0={}
        self.template_1={}
        self.colourlaw=[]

        for phase in numpy.unique(template_0['phase']):
            selection=(template_0['phase']==phase)
            self.template_0['%s' % phase]=template_0[selection]['wave','flux']

        for phase in numpy.unique(template_1['phase']):
            selection=(template_1['phase']==phase)
            self.template_1['%s' % phase]=template_1[selection]['wave','flux']

        return


def compareSALTsurfaces(surface):

    import matplotlib.pyplot as plt

    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)
    
    # -----------  Read in the SALT models -------------------

    surface1=readSALTsurface(JLA.get_full_path(params['orig'])+'salt2-4/')
    surface2=readSALTsurface(JLA.get_full_path(params['comp'])+'salt2-4/')
    

    # -----------  Plot the surfaces ----------------------
    fig1=plt.figure()
    for axes,x1 in enumerate([-2,0,2]):
        ax=fig1.add_subplot(3,1,axes+1)
        flux1=surface1.template_0[options.phase]['flux'] + x1 * surface1.template_1[options.phase]['flux']
        flux2=surface2.template_0[options.phase]['flux'] + x1 * surface1.template_1[options.phase]['flux']
        ax.plot(surface1.template_0[options.phase]['wave'],flux1)
        ax.plot(surface2.template_0[options.phase]['wave'],flux2)

    # -----------  Plot the colour laws ----------------------

    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    wave=numpy.arange(2800.,7000.,10.)
    wave_min=2800.
    wave_max=7000.
    
    #wave_min_reduced=reduced_lambda(wave_min)
    #wave_max_reduced=reduced_lambda(wave_max)

    # See saltextinction.h
    reduced_wave=reduced_lambda(wave)
    reduced_wave_B=reduced_lambda(wave_B)
    # Note that reduced_wave_B=0

    alpha=1.0
    # It's not quite correct
    C=-0.1#??

    for coeff in surface1.colour_law['coeff']:
        alpha-=coeff

    p=alpha * reduced_wave
    for exponent,coeff in enumerate(surface1.colour_law['coeff']):
        p+=coeff*reduced_wave**(exponent+2)

    A_wave=p*C
    ax2.plot(wave, A_wave,label='original')

    alpha=1.0
    for coeff in surface2.colour_law['coeff']:
        alpha-=coeff
    p=alpha * reduced_wave
    for exponent,coeff in enumerate(surface2.colour_law['coeff']):
        p+=coeff*reduced_wave**(exponent+2)

    A_wave=p*C
    ax2.plot(wave, A_wave, label='adjusted')

    # Plot CCM R_V=3.1
    E_BV=0.1
    R_V=3.1
    a_wave=E_BV * R_V * CCM(wave, R_V)
    a_B=E_BV * R_V * CCM(numpy.array([wave_B]),R_V)
    ax2.plot(wave,a_wave-a_B,label='CCM R_V=3.1')

    #CCM R_V=1.0
    R_V=1.0
    a_wave=E_BV * R_V * CCM(wave, R_V)
    a_B=E_BV * R_V * CCM(numpy.array([wave_B]),R_V)
    ax2.plot(wave,a_wave-a_B,label='CCM R_V=1.0')

    # F99 R_V=3.1
    a_wave=E_BV * R_V * Fitz99(wave)
    a_B=E_BV * R_V * Fitz99(numpy.array([wave_B]))
    ax2.plot(wave,a_wave-a_B,label='F99 R_V=3.1')
    
    ax2.legend()

    plt.show()
    plt.close()

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="SALTmodels.config",
                      help="Directories containing the SALT models")

    parser.add_option("-p", "--phase", dest="phase", default=0.0,
                      help="Lightcurve phase")

    (options, args) = parser.parse_args()


    compareSALTsurfaces(options)
