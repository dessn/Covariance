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

def modify_CfA4(dst):
    f=open(dst)
    lines=f.readlines()
    f.close()

    # We overwrite the file
    # You may want to do more checks than you do"
    f=open(dst,"w")
    for line in lines:
        if 'KEPLERCAM::V' in line:
            f.write(line.replace("::V","::Vc"))
        elif 'KEPLERCAM::r' in line:
            f.write(line.replace("::r","::rc"))
        elif 'KEPLERCAM::i' in line:
            f.write(line.replace("::i","::ic"))
        else: 
            f.write(line)
    f.close()

    return


surveys=OrderedDict()

surveys['CSP+CfA']={'input':'nearby/CfA-CSP_overlap/concat','output':'CfA-CSP_overlap'}
surveys['CSP']={'input':'nearby/CSP/DS17','output':'CSP'}
surveys['CfA3']={'input':'nearby/CfA3/DS17','output':'CfA3'}
surveys['CfA4']={'input':'nearby/CfA4/DS17','output':'CfA4'}
surveys['DES3yr']={'input':'DES3YR/SMPv5','output':'DES'}

params=JLA.build_dictionary(options.config)

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

    p2=['sn2009hp',
        'sn2009ig',
        'sn2009jr',
        'sn2009kk',
        'sn2009kq',
        'sn2009le',
        'sn2009lf',
        'sn2009li',
        'sn2009na',
        'sn2009nq',
        'sn2010A',
        'sn2010H',
        'sn2010Y',
        'sn2010ag',
        'sn2010ai',
        'sn2010cr',
        'sn2010dt',
        'sn2010dw',
        'snPTF10bjs']


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
            # For CfA4 SNe, modify the the filter names if the SN was affected by the pre-bakeout 
            if survey=="CfA4" and SN in p2:
                modify_CfA4(dst)

file.close()
