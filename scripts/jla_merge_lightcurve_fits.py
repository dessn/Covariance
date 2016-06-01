"""Pythom program to merge the lightcurve fit results into a single format"""

from optparse import OptionParser

def merge_lightcurve_fits(options):
    """Pythom program to merge the lightcurve fit results into a sigle format"""
    import numpy
    import astropy
    import JLA_library as JLA
    from astropy.table import Table, MaskedColumn, vstack

    params = JLA.build_dictionary(options.config)

    # ---------------- JLA ------------------------
    lightCurveFits = JLA.get_full_path(params['JLAlightCurveFits'])
    f=open(lightCurveFits)
    header=f.readlines()
    f.close()
    names=header[0].strip('#').split()

    # I imagine that the tables package in astropy could also be used to read the ascii input file
    SNeSpec = Table(numpy.genfromtxt(lightCurveFits,
                               skip_header=1,
                               dtype='S12,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8',
                               names=names))

    nSNeSpec=len(SNeSpec)
    print 'There are %d SNe from the spectrscopically confirmed sample' % (nSNeSpec)

    # Add an extra column to the table
    SNeSpec['source']=['JLA']*nSNeSpec

    # ---------------- Shuvo's sample ------------------------
    # Photometrically identified SNe in Shuvo's sample, if the parameter exists
    if params['photLightCurveFits']!='None':
        lightCurveFits = JLA.get_full_path(params['photLightCurveFits'])
        SNePhot=Table.read(lightCurveFits, format='fits')
        nSNePhot=len(SNePhot)

        print 'There are %d SNe from the photometric sample' % (nSNePhot)

        # Converting from Shuvo's names to thosed used by JLA
        conversion={'name':'name_adj', 'zcmb':None, 'zhel':'z', 'dz':None, 'mb':'mb', 'dmb':'emb', 'x1':'x1', 'dx1':'ex1', 'color':'c', 'dcolor':'ec', '3rdvar':'col27', 'd3rdvar':'d3rdvar', 'tmax':None, 'dtmax':None, 'cov_m_s':'cov_m_x1', 'cov_m_c':'cov_m_c', 'cov_s_c':'cov_x1_c', 'set':None, 'ra':None, 'dec':None, 'biascor':None}

        # Add the uncertainty in the mass column
        SNePhot['d3rdvar']=(SNePhot['col29']+SNePhot['col28'])/2. - SNePhot['col27']

        # Remove columns that are not listed in conversion
    
        for colname in SNePhot.colnames:
            if colname not in conversion.values():
                SNePhot.remove_column(colname)

    
        for key in conversion.keys():
            # Rename the column if it does not already exist
            if conversion[key]!=None and conversion[key]!=key:
                SNePhot.rename_column(conversion[key], key)
            elif conversion[key]==None:
                # Create it, mask it, and fill all values
                SNePhot[key]=MaskedColumn(numpy.zeros(nSNePhot), numpy.ones(nSNePhot,bool))
                SNePhot[key].fill_value=-99 # does not work as expected, so we set it explicitly in the next line
                SNePhot[key]=-99.9
            else:
                # Do nothing if the column already exists
                pass

        # Add the source column
        SNePhot['source']="Phot_Uddin"       

    # ----------------------  CfA4 ----------------------------------
    if params['CfA4LightCurveFits']!='None':
        lightCurveFits = JLA.get_full_path(params['CfA4LightCurveFits'])
        f=open(lightCurveFits)
        header=f.readlines()
        f.close()
        names=header[0].strip('#').split(',')    

        SNeCfA4=Table(numpy.genfromtxt(lightCurveFits,
                                       skip_header=1,
                                       dtype='S12,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8',
                                       names=names,delimiter=','))

        nSNeCfA4=len(SNeCfA4) 
    
        conversion={'name':'name', 'zcmb':None, 'zhel':'z', 'dz':None, 'mb':'mb', 'dmb':'emb', 'x1':'x1', 'dx1':'ex1', 'color':'c', 'dcolor':'ec', '3rdvar':None, 'd3rdvar':None, 'tmax':None, 'dtmax':None, 'cov_m_s':'cov_m_x1', 'cov_m_c':'cov_m_c', 'cov_s_c':'cov_x1_c', 'set':None, 'ra':None, 'dec':None, 'biascor':None}

        # Remove columns that are not listed in conversion
    
        for colname in SNeCfA4.colnames:
            if colname not in conversion.values():
                SNeCfA4.remove_column(colname)
    
        for key in conversion.keys():
            # Rename the column if it does not already exist
            if conversion[key]!=None and conversion[key]!=key:
                SNeCfA4.rename_column(conversion[key], key)
            elif conversion[key]==None:
                # Create it, mask it, and fill all values
                SNeCfA4[key]=MaskedColumn(numpy.zeros(nSNeCfA4), numpy.ones(nSNeCfA4,bool))
                SNeCfA4[key].fill_value=-99 # does not work as expected, so we set it explicitly in the next line
                SNeCfA4[key]=-99.9
            else:
                # Do nothing if the column already exists
                pass

        # Add the source column
        SNeCfA4['source']="CfA4"   

    try:
        SNe=vstack([SNeSpec,SNePhot,SNeCfA4])
    except:
        SNe=SNeSpec

    # Write out the result as a FITS table
    date = JLA.get_date()
    SNe.write('%s_%s.fits' % (options.output, date), format='fits')

    return

if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA parameters")

    PARSER.add_option("-o", "--output", dest="output", default="JLA++",
                      help="Light curve fit parameters")

    (OPTIONS, ARGS) = PARSER.parse_args()

    merge_lightcurve_fits(OPTIONS)
