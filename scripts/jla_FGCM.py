# Program that propogates the uncertainty in the filter curves to the ZP 
import JLA_library as JLA
import numpy
import pyfits as fits

def prop_unc(params,filt):

    # Use the filterwheel to find the filename of the filter
    filterDir=params['DES_instrument']
    filterWheel=JLA.get_full_path(filterDir)+'/Filterwheel'
    filterNames=numpy.genfromtxt(filterWheel,comments='#',usecols=(0,1),dtype='S1,S30',
                                 names=['filterName','filename'])
    filter_filename=filterNames[filterNames['filterName']==filt['filter'][-1:]]['filename'][0]

    # Read in the filter curve
    filterCurve=JLA.filterCurve(JLA.get_full_path(filterDir)+'/'+filter_filename)

    # Set the magnitude of the filter offset
    offset=filt['wavelength']*10.  

    # We compute a number of integrals. First with the filtercurves as is, then with an offset added to the filter curve
    # i) The I0 integral 
    error_I0=2.5 * numpy.log10(filterCurve.I0(0.0)/filterCurve.I0(offset))
    # ii) The chromatic correction.
    # Assumed to be zero for now
    # If the standard filter transmission curves are shifted by 5nm, then all the filters will be out by that same amount
    # This may mean, that the offset is quite small
    #mean_wavelength=filterCurve.mean()
    #I10_std=filterCurve.I1(mean_wavelength,0.0) / filterCurve.I0(0.0)
    #I10_std_offset=filterCurve.I1(mean_wavelength,10.0) / filterCurve.I0(10.0)

    error_chromatic=0.0

    # iii) The error in the AB offset
    # We use the standard filter curve to compute the AB offset
    calspec=JLA.spectrum(fits.getdata(JLA.get_full_path(params['calspec']),1),'CALSPEC')
    error_AB=filterCurve.AB(calspec)-filterCurve.AB(calspec,offset)
    return error_I0,error_chromatic,error_AB
