"""Python program that modify light curves so that theu use the magnitude system developed for S18
"""

# 

# Usage
# jla_modify_lightcurves_S18.py -c configFile
#

from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import JLA_library as JLA
import os
import shutil as sh

# CfA1
# CfA2

modifications={'Riess1999_LC':'CFA1_S18','Jha2006_LC':'CFA2_S18'}
#modifications={}

def modify(lc,options,params):
    f_input=open('%s/%s' % (params['inputDir'],lc))
    lines=f_input.readlines()
    f_input.close()

    mod=False
    
    f_output=open('%s/%s' % (params['outputDir'],lc),'w')
    for line in lines:
        entries=line.split()
        if line[0] in '@#' and len(line)>0:
            # Modify the extinction if required
            if '@MWEBV'== entries[0]:
                if options.adjustExtinction:
                    f_output.write('@MWEBV %6.5f\n' % (float(entries[1])*0.86))
                else:
                    f_output.write(line)
            
            elif '@SURVEY'== entries[0]:
                f_output.write(line)
                if entries[1] in modifications.keys() and not options.sameSystem:
                    mod=True
                    survey=entries[1]

            else:
                f_output.write(line)
            
        else:
            if mod and len(line.strip())>0:
                if survey in modifications.keys():
                    f_output.write('%s %s\n' % (' '.join(entries[:-1]),modifications[survey]))
            else:
                f_output.write(line)
                

    if mod:
        print("Modified %s" % lc)
        
    f_output.close()
    
    return None

def modify_lightcurves_S18(options):
    # Read in the configuration file
    # The configuraiton file contains the location of various files
    params=JLA.build_dictionary(options.config)

    # Make the output directory
    try:
        os.mkdir(params['outputDir'])
    except:
        pass
    
    f=open("%s/%s" % (params['inputDir'],params['inputList']))
    listing=f.readlines()
    f.close()

    for entry in listing:
        items=entry.split()
        if items[1]!='LC':
            sh.copy('%s/%s' % (params['inputDir'],items[2]), '%s/%s' % (params['outputDir'],items[2]))
        else:
            modify(items[2],options,params)

    # Copy accross the covariance files
    for entry in os.listdir(params['inputDir']):
        if 'covmat' in entry:
            sh.copy('%s/%s' % (params['inputDir'],entry), '%s/%s' % (params['outputDir'],entry))

    # Copy accross the listing
    sh.copy('%s/%s' % (params['inputDir'],params['inputList']), '%s/%s' % (params['outputDir'],params['inputList']))
    
            
    return    


if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="S18.config",
                      help="Parameter file containting the location of various JLA files")

    parser.add_option("-e", "--adjustExtinction", dest="adjustExtinction", default=False, action='store_true',
                      help="Adjust the extinction")

    # For backwards compatibility
    parser.add_option("-s", "--sameSystem", dest="sameSystem", default=False, action='store_true',
                      help="Same instrumental files")


    (options, args) = parser.parse_args()


    modify_lightcurves_S18(options)
