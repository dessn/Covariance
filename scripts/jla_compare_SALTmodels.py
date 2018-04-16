"""Python program to compute the difference between SALT models
"""

from optparse import OptionParser
import os
import JLA_library as JLA
import numpy as np

# Usage
# JLA_computeCcal.py optioons
# 
# The method follows the description in section 5.4 of B14

def nmad(data):
    return 1.4826 * np.nanmedian(np.abs(data-np.nanmedian(data)))

def runSALT(SALTpath, SALTmodel, salt_prefix, inputFile, SN):
    
    # Set up the path to the SALT model and the name of the outputFile
    os.environ['SALTPATH']=JLA.get_full_path(SALTpath)
    outputFile=JLA.get_full_path(options.workArea)+'/'+SN+'/'+SN+'_'+SALTmodel+'.dat'
    if os.path.isfile(outputFile):
        pass
    else:
        # Otherwise, do the fit with the date of Max set to the value in the lightcurve file
        print 'Fitting %s using  %s' % (SN,os.environ['SALTPATH'])
        JLA.fitLC(inputFile, outputFile, salt_prefix, forceDayMax=False)
    return outputFile

def compareSALTmodels(options):

    import numpy
    import astropy.io.fits as fits
    from astropy.table import Table
    import matplotlib.pyplot as plt

    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)

    # ---------- Read in the SNe lists -------------------------

    SNeList1 = Table(numpy.genfromtxt(params['data1'],
                                     usecols=(0, 2),
                                     dtype='S30,S100',
                                     names=['id', 'lc']))

    SNeList2 = Table(numpy.genfromtxt(params['data2'],
                                     usecols=(0, 2),
                                     dtype='S30,S100',
                                     names=['id', 'lc']))

    for i,SN in enumerate(SNeList1):
        SNeList1['id'][i]=SNeList1['id'][i].replace('lc-', '').replace('.list', '').replace('_smp', '')

    for i,SN in enumerate(SNeList2):
        SNeList2['id'][i]=SNeList2['id'][i].replace('lc-', '').replace('.list', '').replace('_smp', '')

    salt_prefix=''

    C=[]
    dM=[]
    dX=[]
    dC=[]
    z=[]

    try:
        os.mkdir(options.workArea)
    except:
        pass

    for SN1,SN2 in zip(SNeList1,SNeList2):
        try:
            os.mkdir(options.workArea+'/'+SN1['id'])
        except:
            pass
        
        # Fit the SNe with one model, retrieve the results
        SALTmodel='model1'
        model1=runSALT(params[SALTmodel],SALTmodel,salt_prefix,SN1['lc'],SN1['id'])

        # Fit the SNe with the second model, retrieve the results
        SALTmodel='model2'
        model2=runSALT(params[SALTmodel],SALTmodel,salt_prefix,SN2['lc'],SN2['id'])

        result=JLA.computeOffsets(model1,model2,getRedshift=True)
        dM.append(result[0])
        dX.append(result[1])
        dC.append(result[2])
        z.append(result[3])
        print '%20s at z=%4.3f dM: %7.4f\tdX: %7.4f\tdC:%7.4f' % (SN1['id'],z[-1],dM[-1],dX[-1],dC[-1])
        #C.extend([dM,dX,dC])
        
    diff=Table([dM,dX,dC,z],names=['dM','dX','dC','z'])
    date=JLA.get_date()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    alpha=0.141
    beta=3.101
    dmu=diff['dM']+alpha*diff['dX']-beta*diff['dC']
    ax.plot(diff['z'],dmu,'ro')
    ax.set_xlabel('redshift')
    ax.set_ylabel('$\Delta \mu$')
    ax.set_title(options.output)
    ax.set_ylim(-0.3,0.3)
    print 'NMAD %6.3f off z<0.1 %6.3f off z>0.1 %6.3f' % (nmad(dmu),np.median(dmu[diff['z']<0.1]),np.median(dmu[diff['z']>0.1]))

    plt.savefig('%s.png' % (options.output))
    plt.show()
    plt.close()

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="SALTmodels.config",
                      help="Dirctories containing the SALT models")

    parser.add_option("-w", "--workArea", dest="workArea", default="workArea",
                      help="Work area that stores the light curve fits")

    parser.add_option("-o", "--output", dest="output", default="diff",
                      help="Output file")

    (options, args) = parser.parse_args()


    compareSALTmodels(options)
