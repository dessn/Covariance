import subprocess as sp
import numpy
import os
from scipy import interpolate
import time

class filterCurve:
    """A filter"""
    def __init__(self,arg1,arg2=None):
        if arg2==None:
            f=open(arg1,'r')
            wave=[]
            trans=[]
            for line in f.readlines():
                if line[0]!='#':
                    entries=line.split()
                    wave.append(float(entries[0])) # We use Angstroms for the wavelength in the filter transmission file
                    trans.append(float(entries[1]))
            f.close()
            self.wave=numpy.array(wave)        
            self.trans=numpy.array(trans)
        else:
            self.wave=numpy.array(arg1)        
            self.trans=numpy.array(arg2)

        return

    def flux(self,spectrum,offset=0):
        step=1.0
        wave=numpy.arange(numpy.min(self.wave),numpy.max(self.wave),step)
        trans=interpolate.interp1d(self.wave,self.trans)(wave)
        flux=interpolate.interp1d(spectrum.wave,spectrum.flux)(wave+offset*10.)
        return (trans*flux*wave).sum()  # The extra lambda term arrises because we are working with photon counting devices

    def I0(self,offset=0):
        # Computes the I0 integral in Burke et al.
        step=1.0
        wave=numpy.arange(numpy.min(self.wave),numpy.max(self.wave),step)+offset
        trans=interpolate.interp1d(self.wave,self.trans)(wave-offset)
        return (trans / wave).sum()

    def I1(self,mean,offset=0):
        # Computes the I1 integral in Burke et al.
        step=1.0
        wave=numpy.arange(numpy.min(self.wave),numpy.max(self.wave),step)+offset
        trans=interpolate.interp1d(self.wave,self.trans)(wave-offset)
        return (trans * (wave-mean) / wave).sum()

    def mean(self):
        # Computes the mean photon wavelength
        # Eq. 12 in Burke et al.
        step=1.0
        wave=numpy.arange(numpy.min(self.wave),numpy.max(self.wave),step)
        trans=interpolate.interp1d(self.wave,self.trans)(wave)
        return trans.sum() / (trans / wave).sum()

    def AB(self,spectrum,offset=0):
        # Offset is a wavelength offset in the filter
        step=1.0
        wave=numpy.arange(numpy.min(self.wave),numpy.max(self.wave),step)+offset
        trans=interpolate.interp1d(self.wave,self.trans)(wave-offset)
        flux=interpolate.interp1d(spectrum.wave,spectrum.flux)(wave)
        # AB = 0, corresponds to a flux density of 3631 Jy
        # 1 Jy = 1e-23 erg/s/cm^2/Hz
        # We set Plancks constant to 1 as this factors out
        # h=6.62607e-27 # erg.s
        h=1.0
        c=2.99792e18  # Angstrom/s
        # With the lux of the objects in f_lambda
        int1=flux * trans * (wave) / h / c
        int2= 3631e-23 * trans / (wave) / h
        return -2.5*numpy.log10(int1.sum() / int2.sum())

    def eff(self):
        # Assume a flat spectrum
        step=1.0
        wave=numpy.arange(numpy.min(self.wave),numpy.max(self.wave),step)
        trans=interpolate.interp1d(self.wave,self.trans)(wave)
        int1=trans * wave
        int2=trans
        return int1.sum() / int2.sum()

class spectrum:
    """A spectrum"""
    def __init__(self,spectrum,format='ASCII'):
        wave=[]
        flux=[]
        if format=='ASCII':
            file=open(spectrum,'r')
            for line in file.readlines():
                if line[0]!='#':
                    entries=line.split()
                    wave.append(float(entries[0]))
                    flux.append(float(entries[1]))
            file.close()
            self.wave=numpy.array(wave)        
            self.flux=numpy.array(flux)
        else:
            self.wave=spectrum['WAVELENGTH']        
            self.flux=spectrum['FLUX']
        return


class components:
    def __init__(self,filename):
        f=open(filename)
        data=f.readlines()
        f.close
        self.basis=[]
        self.nPhases=[]
        self.nWaveBins=[]
        self.startPhase=[]
        self.endPhase=[]
        self.startWave=[]
        self.endWave=[]
        self.data=[]

        # For each phase, there are 100 points
        # We store the pca as a numpy array that is indexed as phase,wavelength

        for record in data:
            if record.strip()!='':
                entries=record.split()
                if entries[0]=='BSplineAdaptedBasis':
                    print 'There are %d entries' % len(entries)
                    self.basis.append(entries[0])
                    self.nPhases.append(int(entries[1]))
                    self.nWaveBins.append(int(entries[2]))
                    self.startPhase.append(float(entries[3]))
                    self.endPhase.append(float(entries[4]))
                    self.startWave.append(float(entries[5]))
                    self.endWave.append(float(entries[6]))
                    self.data.append(numpy.array([float(entry) for entry in entries[7:]]).reshape(self.nPhases[-1],self.nWaveBins[-1],order="F"))

                    
def reindex_SNe(snlist, data):
    """ list of indices to reindex the data so that it matches the list of SNe """
    indices = []
    for sn in snlist:
        indices+=[i for i in range(len(data)) if data['name'][i].strip() == sn.strip()]
    return indices
                    
def get_date():
    ymd = time.gmtime(time.time())
    return '%4d%02d%02d' % (ymd[0], ymd[1], ymd[2])

def get_full_path(path):
    """If required, return to the full path"""

    if '$' in path:
        rootDir=os.environ[path[path.index('$')+1:path.index('/')]]
        return rootDir+path[path.index('/'):]
    else:
        return path
    
def survey(sn):
    """ assigns SN to sample; works for both formats of sn['name'] (list or str)
    """
    # This will eventually fail as SN names will change with the survey
    if len(sn['name']) > 1:
        name = sn['name']
    else:
        name = sn['name'][0]
    SNLS=['D1', 'D2', 'D3', 'D4']
    if name[0:3] == 'DES':
        return 'DES'
    if name[0:4] == 'SDSS':
        return 'SDSS'
    elif name[2:4] in SNLS or name[0:] in SNLS:
        return 'SNLS'
    elif name[0:2] == 'sn':
        return 'nearby'
    else:
        return 'high-z'

def build_dictionary(f):
    """
    Creates a dictionary from an input file
    Input file contains two columns, the parameter and the parameter value
    If required, environment variables are 
    """
    f=open(f,'r')
    lines=f.readlines()
    f.close()
    parameters={}

    for line in lines:
        if len(line.strip()) > 0: # Not a blank line
            if line[0]!='#':
                entries=line.split()
                parameters[entries[0]]=entries[1]

    return parameters

def SALTmodels(saltModels):
    
    f=open(saltModels)
    lines=f.readlines()
    f.close()
    
    d=[]
    o=[]
    for line in lines:
        entries=line.split('\t')
        d.append(entries[0].strip())
        o.append('')

    # We use a structured array to record the SALT models
    models=numpy.zeros(len(d),dtype=[('offset','a10'),('directory','a10')])

    models['offset']=o
    models['directory']=d

    return models

def fitLC(inputFile, outputFile, salt_prefix='',forceDayMax=False, salt_path=None):
    """fits the lightcurve with the specified SALT2.4 lightcurve fitter;
    calls snfit (can specify path to executable)
    """

    # salt_path is where the SALT2 model and instrument files are located
    # salt_prefix is the where the snls code is located

    if salt_path!=None:
        os.environ['SALTPATH']=salt_path
    if forceDayMax:
        # Get the data of maximum light
        f=open(inputFile)
        lines=f.readlines()
        f.close()
        for line in lines:
            entries=line.split()
            if entries[0]=="@DayMax":
                DayMax=entries[1]
                break
        cmd = salt_prefix + 'snfit '+inputFile+' -o '+outputFile + " -f DayMax " + DayMax
    else:
        cmd = salt_prefix + 'snfit '+inputFile+' -o '+outputFile
    # One should write any errors to a log file
    FNULL = open(os.devnull, 'w')
    sp.call(cmd, shell=True, stdout=FNULL, stderr=sp.STDOUT)
    FNULL.close()
    return

def modelLC(inputFile, resultFile, modelFile, salt_prefix=''):
    """SALT2 model for lightcurve; calls snlc (can specify path to executable)
    """
    cmd = salt_prefix + 'snlc '+ inputFile+' -p '+resultFile+' -o '+modelFile
    # One should write any errors to a log file
    FNULL = open(os.devnull, 'w')
    sp.call(cmd,shell=True, stdout=FNULL, stderr=sp.STDOUT)
    return

def computeOffsets(nominalResult,perturbedResult,getRedshift=False):
    """ returns offsets in m_B, X_1 and C when the fit is perturbed 
    """
    # Nominal result
    f=open(nominalResult)
    lines=f.readlines()
    f.close()
    for line in lines:
        entries=line.split()
        if entries[0]=='Color':
            C_nominal=float(entries[1])
            break
    for line in lines:
        entries=line.split()
        if entries[0]=='X1':
            X_nominal=float(entries[1])
            break
    for line in lines:
        entries=line.split()
        if line!='\n':
            if entries[0]=='RestFrameMag_0_B':
                M_nominal=float(entries[1])
                break

    # Perturbed result
    f=open(perturbedResult)
    lines=f.readlines()
    f.close()
    for line in lines:
        entries=line.split()
        if entries[0]=='Color':
            C_perturbed=float(entries[1])
            break
    for line in lines:
        entries=line.split()
        if entries[0]=='X1':
            X_perturbed=float(entries[1])
            break
    for line in lines:
        entries=line.split()
        if line!='\n':
            if entries[0]=='RestFrameMag_0_B':
                M_perturbed=float(entries[1])
                break

    if getRedshift:
        for line in lines:
            entries=line.split()
            if line!='\n':
                if entries[0]=='Redshift':
                    z=float(entries[1])
                    break
        return M_perturbed-M_nominal,X_perturbed-X_nominal,C_perturbed-C_nominal,z

    return M_perturbed-M_nominal,X_perturbed-X_nominal,C_perturbed-C_nominal



def getCovMatfile(inputFile):
    f=open(inputFile)
    lines=f.readlines()
    f.close()
    for line in lines:
        entries=line.split()
        if len(entries) > 0:
            if entries[0]=='@COVMAT':
                return entries[1]
            
    return None

def getDateOfMax(inputFile):
    print inputFile
    f=open(inputFile)
    lines=f.readlines()
    f.close()
    for line in lines:
        entries=line.split()
        if len(entries) > 0:
            if entries[0]=='DayMax':
                return entries[1],entries[2]
            
    return None,None



def insertDateOfMax(SN, inputFile, outputFile, force=False, params={}):
    # If the date of Max is already included skip

    try:
        os.mkdir('workArea')
    except:
        pass
    cwd=os.getcwd()
    SNname=SN.replace("lc-","").replace("_smp","").replace(".list","")
    workArea='workArea/'+SNname

    try:
        os.mkdir(workArea)
    except:
        pass
    os.chdir(workArea)

    f=open(inputFile)
    lines=f.readlines()
    f.close()
    dayMax=False
    covMat=None
    for line in lines:
        if "DayMax" in line:
            dayMax=True
        if "COVMAT" in line:
            covMat=line.split()[1]


    if dayMax and not force:
        cmd='cp %s %s' % (inputFile,outputFile)
        print 'Copying the lightcurve for %s' % SNname 
        sp.call(cmd,shell=True)
    else:
        print 'Fitting and adjusting %s' % SNname 
        outputFile1=SNname+'.dat3'
        fitLC(inputFile, outputFile1, salt_prefix='', forceDayMax=False, salt_path=get_full_path(params['defsaltModel']))
        date,error=getDateOfMax(outputFile1)
        # Remove the old date of max and insert the new one
        lc=open(inputFile)
        lc_lines=lc.readlines()
        lc.close()
        lc=open(outputFile,'w')
        lc.write('@DayMax %s %s\n' % (date, error))
        for line in lc_lines:
            if 'DayMax' in line:
                pass
            else:
                lc.write(line)
        lc.close()


    if covMat!=None:
        inputFile2=os.path.split(inputFile)[0]+'/'+covMat
        outputFile2=os.path.split(outputFile)[0]+'/'+covMat
        cmd='cp %s %s' % (inputFile2,outputFile2)
        print 'Copying the covariance matrix for %s' % SN 
        sp.call(cmd,shell=True)

    os.chdir(cwd)
    

    return None

def compute_extinction_offset(SN, inputFile, offset, workArea, salt_path=''):
    # salt_path is where the SALT2 model and instrument files are located

    # Write the result and the covariance matrix to the area you want to work in
    if salt_path!='':
        os.environ['SALTPATH']=salt_path
    print os.environ['SALTPATH']

    salt_prefix=''

    try:
        os.mkdir(workArea)
    except:
        pass
    cwd=os.getcwd()
    workArea=workArea + '/' + SN
    try:
        os.mkdir(workArea)
    except:
        pass
    os.chdir(workArea)
    cmd='cp '+inputFile+' .'
    sp.call(cmd,shell=True)
    covMat=getCovMatfile(inputFile)
    if covMat!=None:
        cmd2='cp '+os.path.split(inputFile)[0]+'/'+covMat+' .'
        try:
            sp.call(cmd2,shell=True)
        except:
            print 'No covariance matrix'


    # Run snfit twice, once with the correct value of E(B-V) and again with the adjusted value
    # We force the fit to use the date of Max that is set in the light curve file

    inputFile1=os.path.split(inputFile)[1]
    outputFile1=SN+'.fit'
    if os.path.isfile(outputFile1):
        print 'Skipping %s, output file %s already exists' % (SN,outputFile1)
    else:
        fitLC(inputFile1, outputFile1, salt_prefix, True, salt_path)

    # Adjust the value of E(B-V)
    inputFile2=inputFile1+'2'
    outputFile2=outputFile1+'2'
    cmd='awk \'{if ($1=="@MWEBV") print $1,$2*%3.1f; else print $0}\' %s > %s' % (1.0+offset,inputFile1,inputFile2)
    sp.call(cmd,shell=True)

    if os.path.isfile(outputFile2):
        print 'Skipping %s, output file %s already exists' % (SN,outputFile2)
    else:
        fitLC(inputFile2, outputFile2, salt_prefix, True, salt_path)

    f=open(outputFile1)
    lines=f.readlines()
    f.close()
    for line in lines:
        entries=line.split()
        if entries[0]=='Color':
            C_nominal=float(entries[1])
            break
    for line in lines:
        entries=line.split()
        if entries[0]=='X1':
            X_nominal=float(entries[1])
            break
    for line in lines:
        entries=line.split()
        if line!='\n':
            if entries[0]=='RestFrameMag_0_B':
                M_nominal=float(entries[1])
                break

    # Perturbed result
    f=open(outputFile2)
    lines=f.readlines()
    f.close()
    for line in lines:
        entries=line.split()
        if entries[0]=='Color':
            C_perturbed=float(entries[1])
            break
    for line in lines:
        entries=line.split()
        if entries[0]=='X1':
            X_perturbed=float(entries[1])
            break
    for line in lines:
        entries=line.split()
        if line!='\n':
            if entries[0]=='RestFrameMag_0_B':
                M_perturbed=float(entries[1])
                break

    os.chdir(cwd)
    return M_perturbed-M_nominal,X_perturbed-X_nominal,C_perturbed-C_nominal

def smooth(x,y,width=11):
    # Compute the median in bins of width objects sorted in redshift
    # This smoothing algorithm differs from that in B14
    z=numpy.argsort(x)
    # This excludes SNe that do not fit into a whole bin at the end
    med1=numpy.median(x[z][:(x.size // width) * width].reshape(-1,width),axis=1)
    med2=numpy.median(y[z][:(y.size // width) * width].reshape(-1,width),axis=1)
    # We compute the median of the last few points and add that
    med1=numpy.append(med1,numpy.median(x[z][-1.0 * (x.size % width):]))
    med2=numpy.append(med2,numpy.median(y[z][-1.0 * (y.size % width):]))

    # Now, add the first and last points
    # To ensure that these points do not have too much weight, we simply set them to have the values
    # of the neighbouring point
    med1=numpy.append(x[z[0]],med1)
    med1=numpy.append(med1,x[z[-1]])
    med2=numpy.append(med2[0],med2)
    med2=numpy.append(med2,med2[-1])
   
    # Fit a spline curve to these points
    ## print len(med1), len(med2) ## too short for cubic spline for high-z sample
    tck = interpolate.splrep(med1, med2, s=0)
    # Examine the noise around the best spline fit
    # residual=(y[z]-interpolate.splev(x[z],tck,der=0))**2.
    result=interpolate.splev(x,tck,der=0)
    forPlotting=[med1,med2]
    return forPlotting,result
 
def smoothJLA():
    # A place holder for an algorithm that is closer to JLA
    # In bins of 10 SNe, compute the scatter in that bin
    # Take the median value of the scatter
    # Adjust the smoothness parameter of the spline so that the summed difference between the 
    # spline and the data squared is the number of SN time the scatter squared.
    return None

def A_fn(alpha,beta, n):
    A = numpy.zeros((n,3*n))
    for i in range(n):
        A[i,3*i] += 1
        A[i,3*i+1] += alpha
        A[i,3*i+2] += -beta
    return A
