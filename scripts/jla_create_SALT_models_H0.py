"""Python program to create the SALT models surfaces from DES using the JLA surfaces 
"""

from optparse import OptionParser
import numpy
import shutil
import JLA_library as JLA
from astropy.table import Table

def mkdir(directory):
    import os
    try:
        os.mkdir(directory)
    except:
        print '%s aslready exists' % directory

def offsetFilter(filt,instrument):
    print 'Offsetting %s' % filt
    f=Table.read(filt,format='ascii', comment='#')
    if 'comments' in f.meta.keys():
        f.meta['comments'].append('Offset curve by -1nm')
    else:
        f.meta['comments']=['Offset curve by -1nm']
#    f['col1']-=10. What was this for?
    f.write(filt,format='ascii.fast_no_header')
    return

def offsetZP(magsys,filt,instrument,fitmodel):
    # We loose the @ symbol with the Table object, so we do it by editing the file directly
    print 'Modifying %s %s in %s' % (instrument,filt,magsys)
    f=open(magsys,'r')
    data=f.readlines()
    f.close()

    f=open(magsys,'w')
    f.write('# Offset %s ZP by 0.01\n' % filt)
    for line in data:
        if len(line.strip())==0:
            f.write("\n")
        elif line[0] in ['#','@']:
            f.write(line)
        else:
            entries=line.split()
            if entries[0]==instrument and entries[1]==filt:
                f.write('%s %s %7.3f\n' % (entries[0],entries[1],float(entries[2])+0.01))
            else:
                f.write('%s %s %7.3f\n' % (entries[0],entries[1],float(entries[2])))

    f.close()

    return

def offsetZP2(magSys,offset,filt):
    print "Modifying ZP for %s by %5.3f" % (filt,offset)
    f=open(magSys,'r')
    data=f.readlines()
    f.close()

    instrument=filt.split()[0]
    filt=filt.split()[1]

    f=open(magSys,'w')
    f.write('# Offset %s ZP by %5.3f\n' % (filt,offset))
    for line in data:
        if len(line.strip())==0:
            f.write("\n")
        elif line[0] in ['#','@']:
            f.write(line)
        else:
            entries=line.split()
            if entries[0]==instrument and entries[1]==filt:
                f.write('%s %s %7.3f\n' % (entries[0],entries[1],float(entries[2])+offset))
            else:
                f.write('%s %s %7.3f\n' % (entries[0],entries[1],float(entries[2])))

    f.close()
    return
    
def updateKeplercam(options,params,model):
    try:
        shutil.rmtree(options.output+'/'+model['modelNumber']+'/snfit_data/Instruments/Keplercam')
    except:
        pass
    shutil.copytree(JLA.get_full_path(params['KelplerCam_instrument']),options.output+'/'+model['modelNumber']+'/snfit_data/Instruments/Keplercam')
    
    return

def updateCfA1(options,params,model):

    shutil.copytree(JLA.get_full_path(params['CfA1_instrument']),options.output+'/'+model['modelNumber']+'/snfit_data/Instruments/CfA1')
    
    return

def updateCfA2(options,params,model):

    shutil.copytree(JLA.get_full_path(params['CfA2_instrument']),options.output+'/'+model['modelNumber']+'/snfit_data/Instruments/CfA2')
    
    return


def updateDES(options,params,model):
    try:
        shutil.rmtree(options.output+'/'+model['modelNumber']+'/snfit_data/Instruments/DECam')
    except:
        pass
    shutil.copytree(JLA.get_full_path(params['DES_instrument']),options.output+'/'+model['modelNumber']+'/snfit_data/Instruments/DECam')
    
    # Update the DES magnitude system
    shutil.copy(JLA.get_full_path(params['DES_magsys']),options.output+'/'+model['modelNumber']+'/snfit_data/MagSys/')

    return

def create_Models(options):
    import os

    params=JLA.build_dictionary(options.config)

    try:
        os.mkdir(options.output)
    except:
        print "Directory %s already exists" % (options.output)

    # Read in the SALT models that will be kept
    print options.modelList
    SALTmodels=Table.read(options.modelList,format='ascii',names=['ID','Description'],data_start=0)

    # Read in the models for which the magnitude will be adjusted
    try:
        magOffsets=Table.read(options.magOffsetList,format='ascii',names=['Model','Filter','Offset','MagSys'],data_start=1,delimiter='\t',comment='#')
    except:
        magOffsets=[]

    modelList=[]
    for model in os.listdir(JLA.get_full_path(options.base)):
        if model in SALTmodels['ID']:
            print "Copying across %s" % model
            modelList.append(model)
            shutil.copytree(options.base+'/'+model,options.output+'/'+model)
            # Copy salt2 directory to salt2-4
            shutil.copytree(options.output+'/'+model+'/snfit_data/salt2',options.output+'/'+model+'/snfit_data/salt2-4')
            # Update fitmodel.card
            shutil.copy(JLA.get_full_path(params['fitmodel']),options.output+'/'+model+'/snfit_data/fitmodel.card')
            # Update the Keplercam instrument files
            shutil.rmtree(options.output+'/'+model+'/snfit_data/Instruments/Keplercam')
            shutil.copytree(JLA.get_full_path(params['KelplerCam_instrument']),options.output+'/'+model+'/snfit_data/Instruments/Keplercam')
            # Create the CfA1 instrument files
            shutil.copytree(JLA.get_full_path(params['CfA1_instrument']),options.output+'/'+model+'/snfit_data/Instruments/CfA1')
            # Create the CfA2 instrument files
            shutil.copytree(JLA.get_full_path(params['CfA2_instrument']),options.output+'/'+model+'/snfit_data/Instruments/CfA2')
            # Update the Magnitude system
            shutil.copy(JLA.get_full_path(params['Magsys']),options.output+'/'+model+'/snfit_data/MagSys/')
            # This is not needed for CSP as the instrument files and magnitude system have not chnaged since JLA
        else:
            print "Excluding %s" % model


    print 'We start with %d models from JLA' % (len(modelList))

    # ---------  Add new models --------------

    newModels=Table.read(options.add,format='ascii', comment='#')
    for model in newModels:
        # Copy accross the base model
        shutil.copytree(JLA.get_full_path(model['baseModel']),options.output+'/'+model['modelNumber'])
        print 'Creating %s' % (model['modelNumber'])

        # Copy salt2 directory to salt2-4
        shutil.copytree(options.output+'/'+model['modelNumber']+'/snfit_data/salt2',options.output+'/'+model['modelNumber']+'/snfit_data/salt2-4')

        # Remove the old base instrument, if it exists and replace it with a new one
        try:
            shutil.rmtree(options.output+'/'+model['modelNumber']+'/snfit_data/'+model['fitmodel'])
        except:
            pass

        shutil.copytree(JLA.get_full_path(model['baseInstrument']+model['fitmodel']),options.output+'/'+model['modelNumber']+'/snfit_data/'+model['fitmodel'])

        # Remove the old MagSys directory and replace it with the new one
        shutil.rmtree(options.output+'/'+model['modelNumber']+'/snfit_data/MagSys')
        shutil.copytree(JLA.get_full_path(model['baseInstrument'])+'MagSys',options.output+'/'+model['modelNumber']+'/snfit_data/MagSys')
        
        # Replace the fitmodel.card it with the new one
        shutil.copy(JLA.get_full_path(model['baseInstrument'])+'fitmodel.card',options.output+'/'+model['modelNumber']+'/snfit_data/fitmodel.card')

        # Modify filter curve and ZP
        if model['Type']=='filt':
            offsetFilter(options.output+'/'+model['modelNumber']+'/snfit_data/'+model['fitmodel']+'/'+model['Filter'],model['Instrument'])
        else:
            offsetZP(options.output+'/'+model['modelNumber']+'/snfit_data/MagSys/'+model['MagSys'],model['ShortName'],model['Instrument'],model['fitmodel'])

        # We update the instruments in the newly created surfaces
        # We've added new filters for Keplercam and new instruments for CfA1 and CfA2

        # We should try to generalise this, as this will become very complex as more instruments are added.
        if model['Instrument']=='SWOPE2':
            # Update the Keplercam instrument files
            updateKeplercam(options,params,model)
            updateCfA1(options,params,model)
            updateCfA2(options,params,model)
        elif model['Instrument']=='KEPLERCAM':
            updateCfA1(options,params,model)
            updateCfA2(options,params,model)
            # New SWOPE2 filters are already inlcuded in the instrument files
        elif model['Instrument']=='CfA1':
            updateKeplercam(options,params,model)
            updateCfA2(options,params,model)
        elif model['Instrument']=='CfA2':
            updateKeplercam(options,params,model)
            updateCfA1(options,params,model)
            pass

        modelList.append(model['modelNumber'])

    print 'We now have %d models' % (len(modelList))
    
    # ---- Copy accross the saltModels.list ----

    shutil.copy(options.modelList,options.output+'/saltModels.list')
    

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-b", "--base", dest="base", default=None,
                      help="base model surfaces")

    parser.add_option("-m", "--modelList", dest="modelList", default="saltModels.list",
                      help="list of models to keep")

    parser.add_option("-c", "--config", dest="config", default='DES.config',
                      help="configuration file")

    parser.add_option("-a", "--add", dest="add", default=None,
                      help="List of instruments/filters to add")

    parser.add_option("-o", "--output", dest="output", default=None,
                      help="output model")

    (options, args) = parser.parse_args()


    create_Models(options)
