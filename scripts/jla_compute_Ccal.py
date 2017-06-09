"""Python program to compute the C_cal matrices
"""

from optparse import OptionParser
import os
import JLA_library as JLA
import numpy as np

# Usage
# JLA_computeCcal.py optioons
# 
# The method follows the description in section 5.4 of B14


def runSALT(SALTpath, SALTmodel, salt_prefix, inputFile, SN):
    import os
    
    # Set up the path to the SALT model and the name of the outputFile
    os.environ['SALTPATH']=SALTpath+SALTmodel['directory']+'/snfit_data/'
    outputFile=JLA.get_full_path(options.workArea)+'/'+SN+'/'+SN+'_'+SALTmodel['directory']+'.dat'
    if os.path.isfile(outputFile):
        pass
        #print outputFile, "Skipping, fit with SALT model %s for %s already done" % (SALTmodel['directory'],os.path.split(inputFile)[1])
    else:
        # Otherwise, do the fit with the date of Max set to the value in the lightcurve file
        JLA.fitLC(inputFile, outputFile, salt_prefix) #, forceDayMax=True)
        #print inputFile, outputFile
    return outputFile



def compute_Ccal(options):
    """Python program to compute Ccal
    """

    import numpy
    import astropy.io.fits as fits
    from astropy.table import Table
    import multiprocessing as mp
    import matplotlib.pyplot as plt

    # -----------  Read in the configuration file ------------

    params=JLA.build_dictionary(options.config)
    try:
        salt_prefix = params['saltPrefix']
    except KeyError:
        salt_prefix = ''

    # ---------- Read in the SNe list -------------------------

    SNeList = Table(numpy.genfromtxt(options.SNlist,
                                     usecols=(0, 2),
                                     dtype='S30,S100',
                                     names=['id', 'lc']))


    for i,SN in enumerate(SNeList):
        SNeList['id'][i]=SNeList['id'][i].replace('lc-', '').replace('_smp.list', '')

    # ----------  Read in the SN light curve fits ------------
    # This is used to get the SN redshifts which are used in smoothing the Jacbian

    lcfile = JLA.get_full_path(params[options.lcfits])
    SNe = Table.read(lcfile, format='fits')

    # Make sure that the order is correct
    indices = JLA.reindex_SNe(SNeList['id'], SNe)
    SNe = SNe[indices]


    # -----------  Set up the structures to handle the different salt models -------
    # The first model is the unperturbed salt model
    SALTpath=JLA.get_full_path(params['saltPath'])

    SALTmodels=JLA.SALTmodels(SALTpath+'/saltModels.list')
    nSALTmodels=len(SALTmodels)-1
    print SALTmodels, nSALTmodels

    nSNe=len(SNeList)
    print 'There are %d SNe in the sample' % (nSNe)
    print 'There are %d SALT models' % (nSALTmodels)

    # Add a survey column, which we use with the smoothing, and the redshift
    SNeList['survey'] = numpy.zeros(nSNe,'a10')
    SNeList['z'] = SNe['zhel']

    # Identify the SNLS, SDSS, HST and low-z SNe. We use this when smoothing the Jacobian
    # There is rather inelegant 
    # We still need to allow for Vanina's naming convention when doing this for the photometric sample
    for i,SN in enumerate(SNeList):
        if SN['id'][0:4]=='SDSS':
            SNeList['survey'][i]='SDSS'
        elif SN['id'][2:4] in ['D1','D2','D3','D4']:
            SNeList['survey'][i]='SNLS'
        elif SN['id'][0:3] in ['DES']:
            SNeList['survey'][i]='DES'
        elif SN['id'][0:2]=='sn':
            SNeList['survey'][i]='nearby'
        else:
            SNeList['survey'][i]='high-z'

    # -----------   Read in the calibration matrix -----------------
    Cal=fits.getdata(JLA.get_full_path(params['C_kappa']))

    # Multiply the ZP submatrix by 100^2, and the two ZP-offset submatrices by 100,
    # because the magnitude offsets are 0.01 mag and the units of the covariance matrix are mag
    size=Cal.shape[0] / 2
    Cal[0:size,0:size]=Cal[0:size,0:size]*10000.
    Cal[0:size,size:]*=Cal[0:size,size:]*100.
    Cal[size:,0:size]=Cal[size:,0:size]*100.


    # ------------- Create an area to work in -----------------------
    workArea = JLA.get_full_path(options.workArea)
    try:
        os.mkdir(workArea)
    except:
        pass

    # -----------   The lightcurve fitting --------------------------

    firstSN=True
    
    log=open('log.txt','w')

    for i,SN in enumerate(SNeList):
        J=[]
        try:
            os.mkdir(workArea+'/'+SN['id'])
        except:
            pass

        #firstModel=True
        print 'Examining SN #%d %s' % (i+1,SN['id'])

        # Set up the number of processes
        pool = mp.Pool(processes=int(options.processes))
        # runSALT is the program that does the lightcurve fitting
        results = [pool.apply(runSALT, args=(SALTpath,
                                             SALTmodel,
                                             salt_prefix,
                                             SN['lc'],
                                             SN['id'])) for SALTmodel in SALTmodels]
        for result in results[1:]:
            # The first model is the unperturbed model
            dM,dX,dC=JLA.computeOffsets(results[0],result)
            J.extend([dM,dX,dC])
        pool.close() # This prevents to many open files

        if firstSN:
            J_new=numpy.array(J).reshape(nSALTmodels,3).T
            firstSN=False
        else:
            J_new=numpy.concatenate((J_new,numpy.array(J).reshape(nSALTmodels,3).T),axis=0)

        log.write('%d rows %d columns\n' % (J_new.shape[0],J_new.shape[1]))

    log.close()

    # Compute the new covariance matrix J . Cal . J.T produces a 3 * n_SN by 3 * n_SN matrix
    # J=jacobian

    J_smoothed=numpy.array(J_new)*0.0
    J=J_new

    # We need to concatenate the different samples ...
    
    if options.Plot:
        try:
            os.mkdir('figures')
        except:
            pass               

    nPoints={'SNLS':11,'SDSS':11,'nearby':11,'high-z':11,'DES':11} 
    #sampleList=['nearby','DES']
    sampleList=params['smoothList'].split(',')
    if options.smoothed:
        # We smooth the Jacobian 
        # We roughly follow the method descibed in the footnote of p13 of B14
        for sample in sampleList:
            selection=(SNeList['survey']==sample)
            J_sample=J[numpy.repeat(selection,3)]

            for sys in range(nSALTmodels):
                # We need to convert to a numpy array
                # There is probably a better way
                redshifts=numpy.array([z for z in SNeList[selection]['z']])
                derivatives_mag=J_sample[0::3][:,sys]  # [0::3] = [0,3,6 ...] Every 3rd one
                #print redshifts.shape, derivatives_mag.shape, nPoints[sample]
                forPlotting_mag,res_mag=JLA.smooth(redshifts,derivatives_mag,nPoints[sample])
                derivatives_x1=J_sample[1::3][:,sys]
                forPlotting_x1,res_x1=JLA.smooth(redshifts,derivatives_x1,nPoints[sample])
                derivatives_c=J_sample[2::3][:,sys]
                forPlotting_c,res_c=JLA.smooth(redshifts,derivatives_c,nPoints[sample])

                # We need to insert the new results into the smoothed Jacobian matrix in the correct place
                # The Jacobian ia a 3 * n_SN by nSATLModels matrix
                # The rows are ordered by the mag, stretch and colour of each SNe.
                J_smoothed[numpy.repeat(selection,3),sys]=numpy.concatenate([res_mag,res_x1,res_c]).reshape(3,selection.sum()).ravel('F')

                # If required, make some plots as a way of checking 

                if options.Plot:
                    print 'Creating plot for systematic %d and sample %s' % (sys, sample) 
                    fig=plt.figure()
                    ax1=fig.add_subplot(311)
                    ax2=fig.add_subplot(312)
                    ax3=fig.add_subplot(313)
                    ax1.plot(redshifts,derivatives_mag,'bo')
                    ax1.plot(forPlotting_mag[0],forPlotting_mag[1],'r-')
                    ax1.set_ylabel('mag')
                    ax2.plot(redshifts,derivatives_x1,'bo')
                    ax2.plot(forPlotting_x1[0],forPlotting_x1[1],'r-')
                    ax2.set_ylabel('x1')
                    ax3.plot(redshifts,derivatives_c,'bo')
                    ax3.plot(forPlotting_c[0],forPlotting_c[1],'r-')
                    ax3.set_ylabel('c')
                    ax3.set_xlabel('z')
        
                    plt.savefig('figures/%s_sys_%d.png' % (sample,sys))
                    plt.close()

    date=JLA.get_date()


    fits.writeto('J_%s.fits' % (date) ,J,clobber=True) 
    fits.writeto('J_smoothed_%s.fits' % (date), J_smoothed,clobber=True) 

    # Some matrix arithmatic
    # C_cal is a nSALTmodels by nSALTmodels matrix

    # Read in a smoothed Jacobian specified in the options
    if options.jacobian != None:
        J_smoothed=fits.getdata(options.jacobian)
#    else:
#        # Replace the NaNs in your smoothed Jacobian with zero
#        J_smoothed[numpy.isnan(J_smoothed)]=0

    C=numpy.matrix(J_smoothed)*numpy.matrix(Cal)*numpy.matrix(J_smoothed).T
    if options.output==None:
        fits.writeto('C_cal_%s.fits' % (date), numpy.array(C), clobber=True) 
    else:
        fits.writeto('%s.fits' % (options.output),numpy.array(C),clobber=True)

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA files")

    parser.add_option("-p", "--processes", dest="processes", default="2",
                      help="Number of processes, default 2")

    parser.add_option("-w", "--workArea", dest="workArea", default="../workArea",
                      help="Work area that stores the light curve fits")

    parser.add_option("-o", "--output", dest="output", default=None,
                      help="Output file name")

    parser.add_option("-j", "--jacobian", dest="jacobian", default=None,
                      help="Existing smoothed Jacobian to use")

    parser.add_option("-P", "--Plot", dest="Plot", default=False,
                      action='store_true',
                      help="Plot the Jacobian and the smoothed result")

    parser.add_option("-S", "--smoothed", dest="smoothed", default=True,
                      action='store_false',
                      help="Do not smooth the Jacobian")

    parser.add_option("-s", "--SNlist", dest="SNlist",
                      help="List of SN")

    parser.add_option("-l", "--lcfits", dest="lcfits", default="lightCurveFits",
                      help="Key in config file pointing to lightcurve fit parameters")

    (options, args) = parser.parse_args()


    compute_Ccal(options)
