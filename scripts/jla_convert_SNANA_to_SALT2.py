"""Python program that convert light curves 
"""

# We take the SNANA lightcurves provided by Rick and convert them to the SALT2 format
# It then adds 
# i) Adds noise that seems to correlate with the brightness of the hosts
# ii) Estimates the date of B-band peak

# Usage
# jla_convert_SNANA_to_SALT2.py -c configFile
#
# 


from optparse import OptionParser
import numpy as np
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import JLA_library as JLA
import os

daysOfMax=20.


def valid(s,char):
    if len(s.strip())==0:
        return False
    elif s[0]==char:
        return False
    else:
        return True

def check(obs):
    if obs <=3 or obs >= 1024:
        return True
    else:
        return False

def residuals(p,y,x,e,profile):
    if profile=='polynomial':
        return (y-poly(x,p))/e

def poly(x,p):
    # Should change the basis function
    val=0
    for i in range(len(p)):
        val+=p[i]*x**i
    return val

parameter_names=['SURVEY',
                 'SNID',
                 'IAUC',
                 'TYPE',
                 'FILTERS',
                 'PIXSIZE',
                 'NXPIX',
                 'NYPIX',
                 'FAKE',
                 'RA',
                 'DECL',
                 'MWEBV',
                 'MWEBV_ERR',
                 'PEAKMJD',
                 'REDSHIFT_HELIO',
                 'REDSHIFT_FINAL',
                 'HOSTGAL_OBJID', 
                 'HOSTGAL_PHOTOZ',
                 'HOSTGAL_SPECZ',
                 'HOSTGAL_SNSEP', 
                 'HOSTGAL_MAG',   
                 'HOSTGAL_LOGMASS',   
                 'HOSTGAL_SB_FLUXCAL',
                 'NOBS',
                 'NVAR',
                 'VARLIST']


private_parameter_names=['PRIVATE(DES_snid)',
                         'PRIVATE(DES_cand_type)',
                         'PRIVATE(DES_transient_status)',
                         'PRIVATE(DES_ccdnum)',
                         'PRIVATE(DES_numepochs)',
                         'PRIVATE(DES_numepochs_ml)',
                         'PRIVATE(DES_numepochs_ml_Y0)',
                         'PRIVATE(DES_numepochs_ml_Y1)',
                         'PRIVATE(DES_numepochs_ml_Y2)',
                         'PRIVATE(DES_numepochs_ml_Y3)',
                         'PRIVATE(DES_numepochs_ml_Y4)',
                         'PRIVATE(DES_numepochs_ml_Y5)',
                         'PRIVATE(DES_numepochs_ml_Y6)',
                         'PRIVATE(DES_mjd_trigger)',
                         'PRIVATE(DES_mjd_latest_ml)',
                         'PRIVATE(DES_latest_nite_ml)',
                         'PRIVATE(DES_nobs_offccd)',
                         'PRIVATE(DES_nobs_photerr)',
                         'PRIVATE(DES_hostgal_dlr)',
                         'PRIVATE(DES_hostgal_gradient_g)',
                         'PRIVATE(DES_hostgal_gradient_r)',
                         'PRIVATE(DES_hostgal_gradient_i)',
                         'PRIVATE(DES_hostgal_gradient_z)']

class snanaLightCurve:
    """An snana lightcurve"""
    def __init__(self,filename):
        f=open(filename)
        lines=[s.strip() for s in f.readlines() if valid(s,'#')]
        f.close()

        # Setting up the parameters
        self.parameters={}
        self.private_parameters={}
        for line in lines:
            entries=line.split(':')
            if entries[0] in parameter_names:
                self.parameters[entries[0]]=entries[1].strip()
            if entries[0] in private_parameter_names:
                self.private_parameters[entries[0]]=entries[1].strip()

        # Create a structured array based on self.parameters['VARLIST']
        variables=self.parameters['VARLIST'].split()
        # This is not very robust
        self.data=np.zeros(int(self.parameters['NOBS']),dtype=[(variables[0],float),
                                                               (variables[1],'a1'),
                                                               (variables[2],'a2'),
                                                               (variables[3],float),
                                                               (variables[4],float),
                                                               (variables[5],float),
                                                               (variables[6],float),
                                                               (variables[7],float),
                                                               (variables[8],float),
                                                               (variables[9],int),
                                                               (variables[10],float),])

        print 'Examining %s of type %s' % (self.parameters['SNID'],self.parameters['TYPE'])
        i=0
        for line in lines:
            entries=line.split(':')
            if entries[0]=='OBS':
                values=entries[1].split()
                for j,variable in enumerate(variables): 
                    self.data[i][variable]=values[j]
                i+=1
        return

    def clean(self):
        keep=[]
        for observation in self.data:
            keep.append(check(observation['PHOTFLAG']))
        self.keep=np.array(keep)

        return 

    def write(self,output,format):
        required_variables={'SURVEY':'SURVEY',
                            'SN':'SNID',
                            'RA':'RA',
                            'DECL':'DECL',
                            'Z_HELIO':'REDSHIFT_HELIO',
                            'Z_CMB':'REDSHIFT_FINAL',
                            'MWEBV':'MWEBV',
                            'MWEBV_ERR':'MWEBV_ERR',
                            'DayMax':'PEAKMJD'}

        # Need to determine the field
        # Need to determine the date of max

        analysis_variables={'TYPE':'TYPE'}

        if format=='SALT2.4':
            f=open(output,'w')
            f.write('# ! Translated into SALT2 format using jla_convert_lightcurves.py\n')
            f.write('# !\n')
            f.write('# ! Required variables\n')
            for variable in required_variables.keys():
                f.write('@%s %s\n' % (variable,self.parameters[required_variables[variable]].split()[0]))

            # Add the field name
            f.write('@FIELD %s\n' % (self.data['FIELD'][0]))
            # Add the date of Max
            #f.write('@DayMax %10.3f %4.1f\n' % (self.dateofMax,self.dateofMaxError))

            # We do not include any analysis varaiables for now
            f.write('# !\n')
            f.write('# ! Analysis variables\n')
            f.write('# !\n')
            for variable in analysis_variables.keys():
                f.write('@%s %s\n' % (variable,self.parameters[analysis_variables[variable]].split()[0]))
            # Data
            f.write('# ! --------------------------------------\n')
            f.write('#Date :\n')
            f.write('#Flux :\n')
            f.write('#Fluxerr :\n')
            f.write('#ZP :\n')
            f.write('#Filter :\n')
            f.write('#MagSys :\n')
            f.write('#end :\n')
            # Need to compute a ZP
            for obs in self.data[self.keep]:
                f.write('%-9s %10.4f %10.4f %6.3f DECAM::%1s DES-AB\n' % (obs['MJD'],obs['FLUXCAL'],obs['FLUXCALERR'],obs['ZPFLUX'],obs['BAND']))
            
            f.close()

def convert_lightcurves(options):

    # Read in the configuration file
    # The configuraiton file contains the location of various files
    params=JLA.build_dictionary(options.config)

    # Read in the extra variance
    # This depends on the photometric method. It is lower for SMP 
    ## extraVariance=get_extra_variance(JLA.get_full_path(params['extraVariance']),options)

    snanaDir=JLA.get_full_path(params['snanaLightCurves'])
    saltDir=JLA.get_full_path(params['adjLightCurves'])

    try:
        os.mkdir(saltDir)
    except:
        pass
        
    saltDir=saltDir+'DES/'

    try:
        os.mkdir(saltDir)
    except:
        pass

    for lightcurve in os.listdir(snanaDir):
        if '.dat' in lightcurve:
            # Read in the snana file
            lc=snanaLightCurve(snanaDir+lightcurve)
            lightCurveFile=saltDir+lightcurve.replace('des_real','lc-DES').replace('.dat','.list')
            if lc.parameters['TYPE'].split()[0] in ['1','101']:   # Is a SN Ia or a SN Ia?
                lc.clean()                                # Remove bad photometry
                lc.write(lightCurveFile,options.format)   # Write out the resutlt
    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="DES.config",
                      help="Parameter file containting the location of various JLA files")

    parser.add_option("-p", "--plot", dest="plot", default=False, action='store_true',
                      help="Make plots")

    parser.add_option("-f", "--format", dest="format", default="SALT2.4",
                      help="Output format")


    (options, args) = parser.parse_args()


    convert_lightcurves(options)
