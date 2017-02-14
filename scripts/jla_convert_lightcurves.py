"""Python program that convert light curves 
"""

# Currenly, we take the SNANA lightcurves and modify the SALT2 lightcurves that were provided by Rick
# This program does much the same thing, but also adds an additional source of noise that seems to correlate
# with the brightness of the hosts

from optparse import OptionParser
import numpy as np
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import JLA_library as JLA
import os

daysOfMax=20.

# Usage
# jla_convert_lightcurves.py options
# 

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
                 'SNTYPE',
                 'FILTERS',
                 'PIXSIZE',
                 'NXPIX',
                 'NYPIX',
                 'FAKE',
                 'RA',
                 'DECL',
                 'MWEBV',
                 'MWEBV_ERR',
                 'REDSHIFT_HELIO',
                 'REDSHIFT_FINAL',
                 'HOSTGAL_OBJID', 
                 'HOSTGAL_PHOTOZ',
                 'HOSTGAL_SPECZ',
                 'HOSTGAL_SNSEP', 
                 'HOSTGAL_MAG',   
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
        self.data=np.zeros(int(self.parameters['NOBS']),dtype=[(variables[0],float),
                                                               (variables[1],'a1'),
                                                               (variables[2],'a2'),
                                                               (variables[3],float),
                                                               (variables[4],float),
                                                               (variables[5],int),
                                                               (variables[6],float),
                                                               (variables[7],float),
                                                               (variables[8],float),
                                                               (variables[9],float),
                                                               (variables[10],float),
                                                               (variables[11],float),
                                                               (variables[12],float),
                                                               (variables[13],float),
                                                               (variables[14],int),])

        print 'Examining %s of type %s' % (self.parameters['SNID'],self.parameters['SNTYPE'])
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

    def estimateDateOfMax(self,options):
        # A rough estimate of the date of Max
        # Step 1) Use the maximum flux as the first estimate
        # Step 2) Use a parabolic fit around the date of max

        # Step 1)
        # Select the nearest observer frame filter to the B band
        filters=np.array([('g',475),('r',640),('i',776),('z',925)],dtype=[('filtName','a1'),('wavelength',float)])
        redshift=float(self.parameters['REDSHIFT_FINAL'].split()[0])
        observerFrame=450 * (1+redshift)
        selectedFilter=filters[np.argmin(np.abs(filters['wavelength'] - observerFrame))]['filtName']
        selectedData=self.data[(self.data['FLT']==selectedFilter)]
        index=np.argmax(selectedData['FLUXCAL'])
        dateOfMax=selectedData[index]['MJD']
        
        # Step 2)
        # We fit a polynomial to the data with daysOfMax rest frame days of the peakflux to find an approximate date of maximum light
        selection=(selectedData['MJD'] < dateOfMax+daysOfMax*(1+redshift)) & (selectedData['MJD'] > dateOfMax-daysOfMax*(1+redshift))
        if selection.sum() < 3:
            print 'Expanding the range over the which the coarse fit is done'
            selection=(selectedData['MJD'] < dateOfMax+daysOfMax*(1+redshift)*2.) & (selectedData['MJD'] > dateOfMax-daysOfMax*(1+redshift)*2.)

        mjd=selectedData['MJD'][selection]-dateOfMax # reset the horizontal scale
        flux=selectedData['FLUXCAL'][selection]
        fluxerr=selectedData['FLUXCALERR'][selection]
        guess=(200.0,0.0,-1.0)
        plsq=leastsq(residuals, guess, args=(flux,mjd,fluxerr,'polynomial'), full_output=1)
        # Is a 2nd order polynomial the best choice
        self.dateofMax=-1.0 * plsq[0][1] / plsq[0][2] / 2 + dateOfMax
        self.dateofMaxError=5.0 # This needs to be determined more rigously
        if options.plot:
            fig=plt.figure()
            ax=fig.add_subplot(111)
            ax.errorbar(mjd,flux,yerr=fluxerr,fmt='ro')
            x=np.arange(-20,20,1)
            y=poly(x,plsq[0])
            plt.plot(x,y)
            plt.show()
            plt.close()

        return

    def fitDateOfMax(self,lightCurveFile,params):
        # A full salt2 fit
        outputFile=lightCurveFile.replace('.list','.res')
        os.environ['SALTPATH']=JLA.get_full_path(params['defsaltModel'])
        JLA.fitLC(lightCurveFile, outputFile, salt_prefix='')
        self.dateofMax,self.dateofMaxError=JLA.getDateOfMax(outputFile)
        # Remove the old date of max and insert the new one
        lc=open(lightCurveFile)
        lc_lines=lc.readlines()
        lc.close()
        lc=open(lightCurveFile,'w')
        lc.write('@DayMax %s %s\n' % (self.dateofMax,self.dateofMaxError))
        for line in lc_lines:
            if 'DayMax' in line:
                pass
            else:
                lc.write(line)
        lc.close()
        return

    def addNoise(self, extraVariance):
        filters={'g':0,'r':1,'i':2,'z':3}
        # This part of the code needs to be made more efficient
        # We should compute the spline coefficient when we first read the data
        for index,data in enumerate(self.data):
            # Find the scaleFactor, this depends on field and filter
            for eV in extraVariance.keys():
                ##print eV,data['FLT'],data['FIELD']
                if extraVariance[eV]['filter']==data['FLT'] and data['FIELD'] in extraVariance[eV]['fields']:
                    filterIndex=filters[data['FLT']]
                    # Determine the surfaceBrightess
                    sb=float(self.parameters['HOSTGAL_SB_FLUXCAL'].split()[filterIndex])
                    if sb > 0:
                        surfaceBrightness=27.50-np.log10(sb)
                    else:
                        surfaceBrightness=99
                    if surfaceBrightness < extraVariance[eV]['hostnoise']['mag'][0]:
                        scaleFactor=extraVariance[eV]['hostnoise']['mag'][0]
                    elif  surfaceBrightness > extraVariance[eV]['hostnoise']['mag'][-1]:
                        scaleFactor=extraVariance[eV]['hostnoise']['mag'][-1]
                    else:
                        scaleFactor=interp1d(extraVariance[eV]['hostnoise']['mag'],extraVariance[eV]['hostnoise']['scale'])(surfaceBrightness)
                    self.data[index]['FLUXCALERR']*=scaleFactor

                    break

        return

    def write(self,output,format):
        required_variables={'SURVEY':'SURVEY',
                            'SN':'SNID',
                            'RA':'RA',
                            'DECL':'DECL',
                            'Z_HELIO':'REDSHIFT_HELIO',
                            'Z_CMB':'REDSHIFT_FINAL',
                            'MWEBV':'MWEBV',
                            'MWEBV_ERR':'MWEBV_ERR'}

        # Need to determine the field
        # Need to determine the date of max

        analysis_variables={}

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
            f.write('@DayMax %10.3f %4.1f\n' % (self.dateofMax,self.dateofMaxError))

            # We do not include any analysis varaiables for now
            f.write('# !\n')
            f.write('# ! Analysis variables\n')
            f.write('# !\n')
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
                f.write('%-9s %10.4f %10.4f 27.500 DECAM::%1s DES-AB\n' % (obs['MJD'],obs['FLUXCAL'],obs['FLUXCALERR'],obs['FLT']))
            
            f.close()

    def updateExtinction(self,lightCurveFile,extinction):
        selection=(extinction['filename']==os.path.split(lightCurveFile)[1])
        # Remove the extinction and insert the new one
        if selection.sum()==1:
            lc=open(lightCurveFile)
            lc_lines=lc.readlines()
            lc.close()
            lc=open(lightCurveFile,'w')
            lc.write('@MWEBV %6.4f\n' % (extinction[selection]['extinction'][0]))
            for line in lc_lines:
                if 'MWEBV' in line and 'MWEBV_ERR' not in line:
                    pass
                else:
                    lc.write(line)
            lc.close()
        else:
            print 'Not updating the extinction value'
        return

        
def get_extra_variance(inputFile,options):
    f=open(inputFile)
    lines=[s.strip() for s in f.readlines() if valid(s,'#')]
    f.close()

    extraVariance={}
    libid=0
    mag=[]
    scale=[]

    for line in lines:
        entries=line.split()
        if entries[0]=='LIBID:':
            # A new libid
            if libid>0:
                extraVariance[libid]={'filter':filt,
                                      'fields':fields,
                                      'hostnoise':{'mag':mag,'scale':scale}}
                mag=[]
                scale=[]
            libid+=1
        if entries[0]=='BAND:':
            filt=entries[1]
            fields=entries[3].split('+')
        if entries[0]=='HOSTNOISE:':
            mag.append(float(entries[1]))
            scale.append(float(entries[2]))
            

    extraVariance[libid]={'filter':filt,
                          'fields':fields,
                          'hostnoise':{'mag':mag,'scale':scale}}
    return extraVariance

def get_extinction(fileName,options):

    return np.genfromtxt(fileName,dtype=[('filename','S20'),('extinction',float)])


def convert_lightcurves(options):

    # Read in the configuration file
    params=JLA.build_dictionary(options.config)

    # Read in the extra variance
    extraVariance=get_extra_variance(JLA.get_full_path(params['extraVariance']),options)

    # Read in the extinction values
    extinction=get_extinction(JLA.get_full_path(params['extinction']),options)


    snanaDir=JLA.get_full_path(params['snanaLightCurves'])
    saltDir=JLA.get_full_path(params['adjLightCurves'])+'DES/'
    for lightcurve in os.listdir(snanaDir):
        if '.dat' in lightcurve:
            # Read in the snana file
            lc=snanaLightCurve(snanaDir+lightcurve)
            lightCurveFile=saltDir+lightcurve.replace('des_real','lc-DES').replace('.dat','.list')
            if lc.parameters['SNTYPE'].split()[0] in ['1','101']:   # It is a SN Ia or SN Ia?
                lc.clean()                                # Remove bad photometry
                lc.addNoise(extraVariance)                # Add additional variance to the lightcurve points
                lc.estimateDateOfMax(options)             # Does setting an approximate date of max to the light curve fitting done below.
                # Apply cuts
                # lc.applySNCuts()
                # lc.applySamplingCuts()
                lc.write(lightCurveFile,options.format)   # Write out the resutlt
                lc.fitDateOfMax(lightCurveFile,params)    # Get a more precise estimate of the data of peak brightness
                lc.updateExtinction(lightCurveFile,extinction) # Temporary code
    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA files")

    parser.add_option("-p", "--plot", dest="plot", default=False, action='store_true',
                      help="Make plots")

    parser.add_option("-f", "--format", dest="format", default="SALT2.4",
                      help="Output format")


    (options, args) = parser.parse_args()


    convert_lightcurves(options)
