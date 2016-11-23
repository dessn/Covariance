#!/Users/clidman/anaconda/bin/python
# A program to copy files across from JLA/data to JLA/analysis

# Revised to deal with the merged lightcurve files

import shutil
import os
from optparse import OptionParser
import JLA_library as JLA

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


surveys={'DES':{'input':'DES_3yr_spec/forcePhoto/','output':'DES'},
         'CSP':{'input':'nearby/CSP/PS1s','output':'CSP'},
         'CfA3':{'input':'nearby/CfA3/PS1s','output':'CfA3'},
         'CfA4':{'input':'nearby/CfA4/PS1s','output':'CfA4'},
         'CSP+CfA':{'input':'nearby/CfA-CSP_overlap/concat','output':'CfA-CSP_overlap'}}

params=JLA.build_dictionary(options.config)

SN_list=['sn2009hp']

file=open(options.output,'w')

for survey in surveys.keys():
    try:
        os.mkdir(JLA.get_full_path(params['lightCurves'])+surveys[survey]['output'])
    except:
        pass
    for lightCurve in os.listdir(JLA.get_full_path(params['lightCurvesOrig'])+surveys[survey]['input']):
        if 'list' in lightCurve:
            ori=JLA.get_full_path(params['lightCurvesOrig'])+surveys[survey]['input']+'/'+lightCurve
            dst=JLA.get_full_path(params['lightCurves'])+surveys[survey]['output']+'/'+lightCurve
            SN=lightCurve.replace("lc-","").replace(".list","")
            if surveys[survey]['output'] in ['CfA3','CfA4','CSP'] and  SN in SN_list:
                continue
            SN_list.append(SN)
            file.write('%s LC %s\n'% (lightCurve,dst))
            shutil.copy(ori, dst)
            # For DES SN, modify the magsys from AB to DES-AB
            if survey=="DES":
                modify(dst)

file.close()
