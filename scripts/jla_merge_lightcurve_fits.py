"""Pythom program to merge the lightcurve fit results into a single format"""

from optparse import OptionParser

def getKeywords(lc):
    lcFile=open(lc)
    lines=lcFile.readlines()
    lcFile.close()
    keywords={'HOSTGAL_LOGMASS':-99.9,'e_HOSTGAL_LOGMASS':-99.9}
    for line in lines:
        if line.strip()!='END:':
            entries=line.split()
            if entries[0][-1]==':':
                keywords[entries[0][0:-1]]=entries[1]
            if entries[0]=='HOSTGAL_LOGMASS:':
                keywords['e_HOSTGAL_LOGMASS']=entries[3]
                        
    return keywords

def merge_lightcurve_fits(options):
    """Pythom program to merge the lightcurve fit results into a sigle format"""
    import numpy
    import astropy
    import os
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

    # ---------------- Shuvo's sample aka JLA++  ------------------------
    # Photometrically identified SNe in Shuvo's sample, if the parameter exists
    if params['photLightCurveFits']!='None':
        lightCurveFits = JLA.get_full_path(params['photLightCurveFits'])
        SNePhot=Table.read(lightCurveFits, format='fits')
        nSNePhot=len(SNePhot)

        print 'There are %d SNe from the photometric sample' % (nSNePhot)

        # Converting from Shuvo's names to thosed used by JLA
        conversion={'name':'name_adj', 'zcmb':None, 'zhel':'z', 'dz':None, 'mb':'mb', 'dmb':'emb', 'x1':'x1', 'dx1':'ex1', 'color':'c', 'dcolor':'ec', '3rdvar':'col27', 'd3rdvar':'d3rdvar', 'tmax':None, 'dtmax':None, 'cov_m_s':'cov_m_x1', 'cov_m_c':'cov_m_c', 'cov_s_c':'cov_x1_c', 'set':None, 'ra':'col4', 'dec':'col5', 'biascor':None}

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

        # We also need to gather information on the host mass, the host mass uncertainty and the CMB redshift
        CfA4lightcurves=[]
        CfA4_lcDirectories=params['CfA4MassesAndCMBz'].split(',')
        for lcDir in CfA4_lcDirectories:
            listing=os.listdir(JLA.get_full_path(lcDir))
            for lc in listing:
                CfA4lightcurves.append(JLA.get_full_path(lcDir)+lc)

        for i,SN in enumerate(SNeCfA4):
            for lc in CfA4lightcurves:
                if SN['name'][2:] in lc:
                   keywords=getKeywords(lc)
                   SNeCfA4[i]['zcmb']=keywords['REDSHIFT_CMB']
                   SNeCfA4[i]['3rdvar']=keywords['HOSTGAL_LOGMASS']
                   SNeCfA4[i]['d3rdvar']=keywords['e_HOSTGAL_LOGMASS']
                   SNeCfA4[i]['ra']=keywords['RA']
                   SNeCfA4[i]['dec']=keywords['DECL']
                   if SNeCfA4[i]['3rdvar'] < 0:
                       SNeCfA4[i]['3rdvar']=-99.9
                       SNeCfA4[i]['d3rdvar']=-99.9

    try:
        SNe=vstack([SNeSpec,SNePhot,SNeCfA4])
    except:
        SNe=SNeSpec

    # Write out the result as a FITS table
    date = JLA.get_date()
    SNe.write('%s_%s.fits' % (options.output, date), format='fits', overwrite=True)

    return

if __name__ == '__main__':

    PARSER = OptionParser()

    PARSER.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA parameters")

    PARSER.add_option("-o", "--output", dest="output", default="JLA++",
                      help="Light curve fit parameters")

    (OPTIONS, ARGS) = PARSER.parse_args()

    merge_lightcurve_fits(OPTIONS)
