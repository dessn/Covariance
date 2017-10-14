#!/Users/clidman/anaconda/bin/python
# A program to copy files across from $JLA/SpectroSN/data to JLA/analysis

# Revised to deal with the merged lightcurve files

# 2017-10-13
# Now using the repository $JLA/SpectroSN/data

import shutil
import os
from optparse import OptionParser
import JLA_library as JLA
from collections import OrderedDict

parser = OptionParser()

parser.add_option("-c", "--config", dest="config", default='DES.config',
                  help="configuration file")

parser.add_option("-o", "--output", dest="output", default='DES.list',
                  help="configuration file")

(options, args) = parser.parse_args()

def modify(dst):
    f=open(dst)
    lines=f.readlines()
    f.close()

    # We overwrite the file
    # You may want to do more checks than you do"
    f=open(dst,"w")
    for line in lines:
        f.write(line.replace(" AB"," DES-AB"))
    f.close()

    return


surveys=OrderedDict()

surveys['CSP+CfA']={'input':'nearby/CfA-CSP_overlap/concat','output':'CfA-CSP_overlap'}
surveys['CSP']={'input':'nearby/CSP/PS1s','output':'CSP'}
surveys['CfA3']={'input':'nearby/CfA3/PS1s','output':'CfA3'}
surveys['CfA4']={'input':'nearby/CfA4/PS1s','output':'CfA4'}
surveys['DES3yr']={'input':'DES_3yr_spec/SMPv1','output':'DES'}

params=JLA.build_dictionary(options.config)

#SN_list=['sn2009hp'] Why did we want to skip this SNe?
SN_list=[]

file=open(options.output,'w')

for survey in surveys.keys():
    # Make the subdirectories
    print 'Examining %s' % (survey) 
    try:
        os.mkdir(JLA.get_full_path(params['lightCurves']))
    except:
        pass

    try:
        os.mkdir(JLA.get_full_path(params['adjLightCurves']))
    except:
        pass

    try:
        os.mkdir(JLA.get_full_path(params['lightCurves'])+surveys[survey]['output'])
    except:
        pass

    try:
        os.mkdir(JLA.get_full_path(params['adjLightCurves'])+surveys[survey]['output'])
    except:
        pass

    # Copy the light curves accross and update the list of SNe
    for lightCurve in os.listdir(JLA.get_full_path(params['lightCurvesOrig'])+surveys[survey]['input']):
        if 'list' in lightCurve:
            ori=JLA.get_full_path(params['lightCurvesOrig'])+surveys[survey]['input']+'/'+lightCurve
            dst=JLA.get_full_path(params['lightCurves'])+surveys[survey]['output']+'/'+lightCurve
            SN=lightCurve.replace("lc-","").replace(".list","")
            if surveys[survey]['output'] in ['CfA3','CfA4','CSP'] and  SN in SN_list:
                # If the SN has already been added, then skip
                continue
            SN_list.append(SN)
            file.write('%s LC %s\n'% (lightCurve,dst))
            shutil.copy(ori, dst)
            # For DES SN, modify the magsys from AB to DES-AB
            if survey=="DES":
                modify(dst)

file.close()
