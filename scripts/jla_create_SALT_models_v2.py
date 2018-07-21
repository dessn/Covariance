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
    f['col1']-=10. # Subtract 10nm from filter curve
    f.write(filt,format='ascii.fast_no_header',overwrite=True)
    return

def offsetZP(magsys,filt,instrument,fitmodel):
    # We loose the @ symbol with the Table object, so we do it by editing the file directly
    print 'Modifying %s %s in %s' % (instrument,filt,magsys)
    f=open(magsys,'r')
    data=f.readlines()
    f.close()

    f=open(magsys,'w')
    f.write('# Offset %s %s ZP by 0.01\n' % (instrument,filt))
    for line in data:
        if len(line.strip())==0:
            f.write("\n")
        elif line[0] in ['#','@']:
            f.write(line)
        else:
            entries=line.split()
            if entries[0]==instrument and entries[1]==filt:
                f.write('%s %s %7.3f\n' % (entries[0],entries[1],float(entries[2])+0.01))
                f.write('#%s %s %7.3f\n' % (entries[0],entries[1],float(entries[2])))
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
    shutil.copytree(JLA.get_full_path(params['CfA_instrument']),options.output+'/'+model['modelNumber']+'/snfit_data/Instruments/Keplercam')
    
    # The following is not needed
    ##shutil.copy(JLA.get_full_path(params['CfA_magsys']),options.output+'/'+model['modelNumber']+'/snfit_data/MagSys/')
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
    SALTmodels=Table.read(options.modelList,format='ascii',names=['ID','Type','Instrument','ShortName','fitmodel','MagSys','Filter'],data_start=0)

    modelList=[]
    # Go through the base models
    for model in os.listdir(JLA.get_full_path(options.base)):
        if model in SALTmodels['ID']:

            selection=(SALTmodels['ID']==model)
            modelSel=SALTmodels[selection]

            print "Copying across %s" % model
            modelList.append(model)
            shutil.copytree(options.base+'/'+model,options.output+'/'+model)
            # Copy salt2 directory to salt2-4
            shutil.copytree(options.output+'/'+model+'/snfit_data/salt2',options.output+'/'+model+'/snfit_data/salt2-4')
            # Update fitmodel.card
            shutil.copy(JLA.get_full_path(params['fitmodel']),options.output+'/'+model+'/snfit_data/fitmodel.card')
            # Add the DECam instrument files
            shutil.copytree(JLA.get_full_path(params['DES_instrument']),options.output+'/'+model+'/snfit_data/Instruments/DECam')

            # Update the Keplercam instrument files
            # We added the revised filter curves B,V,r, and i, and Bc, Vc, rc, and ic
            shutil.rmtree(options.output+'/'+model+'/snfit_data/Instruments/Keplercam')
            shutil.copytree(JLA.get_full_path(params['CfA_instrument']),options.output+'/'+model+'/snfit_data/Instruments/Keplercam')

            # Since we overwrite the Keplercam instrument files, we need to offset offset the filter curves
            if modelSel['Type']=='Filter' and modelSel['Instrument']==['KEPLERCAM']:
                print 'Adjusting filter for model %s' % modelSel['ID'][0]
                offsetFilter(options.output+'/'+modelSel['ID'][0]+'/snfit_data/'+modelSel['fitmodel'][0]+'/'+modelSel['Filter'][0],modelSel['Instrument'][0])

            # Add DES magnitude system
            shutil.copy(JLA.get_full_path(params['DES_magsys']),options.output+'/'+model+'/snfit_data/MagSys/')

            # Update the CfA and CSP magnitude systems
            shutil.copy(JLA.get_full_path(params['CfA_magsys']),options.output+'/'+model+'/snfit_data/MagSys/')

            # Since we update the magnitude systems, we need to offset the ZPs for KEPLERCAM, 4SHOOTER, and SWOPE
            if modelSel['Type']=='ZP' and modelSel['Instrument'] in ['KEPLERCAM','4SHOOTER2','SWOPE2']:
                print 'Adjusting ZP for model %s' % modelSel['ID'][0]
                offsetZP(options.output+'/'+modelSel['ID'][0]+'/snfit_data/MagSys/'+modelSel['MagSys'][0],modelSel['ShortName'][0],modelSel['Instrument'][0],modelSel['fitmodel'][0])
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
        print JLA.get_full_path(model['baseInstrument']+model['fitmodel']),options.output+'/'+model['modelNumber']+'/snfit_data/'+model['fitmodel']

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


        # We now update the list of instruments in the newly created surfaces
        # This code is not clear
        if model['Instrument']=='DECAM':
            # Update just the Keplercam instrument files
            # There is no need to update Swope2, as the new filters are in the base model
            updateKeplercam(options,params,model)
        elif model['Instrument']=='KEPLERCAM':
            # Update just the DES instrument files
            # There is no need to update Swope2, as the new filters are in the base model
            updateDES(options,params,model)
        else: # The case for swope ...
            # Update both the Keplercam and DES fles
            updateDES(options,params,model)
            updateKeplercam(options,params,model)




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
