"""
Python program to compute the compute the data of maximum light and add it to the lightcurve file
"""

from optparse import OptionParser

def compute_date_of_max(options):
    import numpy
    from astropy.table import Table
    import JLA_library as JLA

    params=JLA.build_dictionary(options.config)

    # -----------  Read in the configuration file ------------


    lightCurveFits=JLA.get_full_path(params['lightCurveFits'])
    lightCurves=JLA.get_full_path(params['lightCurves'])
    adjlightCurves=JLA.get_full_path(params['adjLightCurves'])
    

    # ---------  Read in the list of SNe ---------------------
    # One can either use an ASCII file with the SN list or a fits file
    if options.SNlist == None:
        SNe = Table.read(lightCurveFits, format='fits')
    else:
    # We use the ascii file, which gives the full path name
        SNe = Table.read(options.SNlist, format='ascii',names=['name','type','lc'],data_start=0)

    nSNe=len(SNe)
    print 'There are %d SNe' % (nSNe)

    # -----------   The lightcurve fitting -------------------

    for SN in SNe:
        if options.SNlist == None:
            SNfile='lc-'+SN['name']+'.list'
            inputFile=lightCurves+SNfile
            outputFile=adjlightCurves+SNfile
        else:
            inputFile=SN['lc']
            outputFile=SN['lc'].replace(lightCurves,adjlightCurves)

        print 'Examining %s' % SN['name']
        print outputFile
        # If needed refit the lightcurve and insert the date of maximum into the input file
        JLA.insertDateOfMax(SN['name'].strip(),inputFile,outputFile,options.force)

    return

# Python program to compute the date of maximum from the lightcurve and to insert it into the lightcurve file

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA parameters")

    parser.add_option("-s", "--SNlist", dest="SNlist", default=None,
                      help="List of SN")

    parser.add_option("-f", "--force", dest="force", default=True, 
                      action='store_true', help="List of SN")

    (options, args) = parser.parse_args()

    compute_date_of_max(options)

