"""Python wrapper for snfit
"""

from optparse import OptionParser
import subprocess as sp
import JLA_library as JLA
import numpy as np
import os
import shutil as sh
import time
c
# Usage
# jla_fitLightCurves.py optioons
# 
# The method follows the procedures specified in e-mail sent by M. Betoule on 15/11/2017

def mkdir(dir):
    try:
        os.mkdir(dir)
    except:
        print("Cannot create %s dir" % dir)

    return

def fitLightCurves(options):

    # Get today's date
    date=date=JLA.get_date()

    # Read in the configuration file
    params=JLA.build_dictionary(options.config)
    
    # Create the directory where the training will be done
    if options.outputDir is not None:
        mkdir(options.outputDir)
    else:
        print("Please specify an output directory")
        exit()
              
    # Read in SN list
    f=open(params['trainingList'])
    listing=f.readlines()
    f.close()

    os.chdir(options.outputDir)
    # Open a log file recording the settings used

    log=open("training.log",'w')
    log.write("# Fitting SN lightcurves \n\n")
    log.write("Date\t%s\n\n" % date)
    log.write("# Parameters\n")
    for key in params.keys():
        log.write("%s\t%s\n" % (key,params[key]))
    log.close()

    if options.fitHostExtinction:
        flags=['-w', '2800', '9000']
    else:
        flags=['-w', '2800', '9000',
               '-f', 'Rv', '3.1',
               '-f', 'Tau', '0.0']
        
    for line in listing:
        entries=line.split()
        SN=entries[0]
        LightCurve=entries[2]
        
        cmd=['snfit',
             LightCurve,
             '-o','%s.dat' % (SN)
             ] + flags

        
        print 'Executing: ',' '.join(cmd)
        sp.call(' '.join(cmd),shell=True)

        
    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="SNFIT.config",
                      help="configuration file containing SALT parameters")

    parser.add_option("-o", "--outputDir", dest="outputDir", default=None,
                      help="Output directory")

    parser.add_option("-f", "--fitHostExtinction", dest="fitHostExtinction", default=False,action="store_true",
                      help="Fit host galaxy extinction")

    

    (options, args) = parser.parse_args()


    fitLightCurves(options)
