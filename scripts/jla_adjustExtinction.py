"""Program to adjust the extinction of the SNe that are used in the training of SALT2
"""

from optparse import OptionParser

def adjustExtinction(options):
    import numpy
    from astropy.table import Table
    import JLA_library as JLA
    import subprocess as sp
    import os

    params=JLA.build_dictionary(options.config)

    # ----------- Correction factor for extinction -----------
    # See ApJ 737 103
    extinctionFactor=0.86

    SNlist=Table.read(params['SNlist'], format='ascii',names=['name','type','lc'],data_start=0)

    try:
        os.mkdir(JLA.get_full_path(params['outputDir']))
    except:
        pass

    for SN in SNlist:
        # Clean the directory
        if options.clean:
            cmd="*fits co* *opt0* *step* re* salt2* spec_coverage* sne_pcafit.list *init* pca_1_opt1_before_final_normalization.list pca_1_superinit.list"
            sp.call(cmd,shell=True)

        # Copy accross the file
        outputDir=JLA.get_full_path(params['outputDir'])
        inputFile=outputDir+os.path.split(SN['lc'])[1]
        cmd='cp %s %s' % (SN['lc'],inputFile)
        sp.call(cmd, shell=True)

        lc=open(inputFile)
        lines=lc.readlines()
        lc.close()

        if SN['type']=='LC':
            # Copy accross the covmat file, if it exists
            covmatExist=False
            for line in lines:
                try:
                    if "@COVMAT"==line.split()[0]:
                        covmat=line.split()[1]
                        covmatExist=True
                        break
                except:
                    pass

            if covmatExist:
                cmd='cp %s %s' % (os.path.split(SN['lc'])[0]+'/'+covmat,outputDir)
                sp.call(cmd, shell=True)
            
            # Adjust the extinction
            mwebv=0.0
            for line in lines:
                if "@MWEBV"==line.split()[0]:
                    mwebv=float(line.split()[1])
                    break

            if mwebv>0:
                # Remove the old extcintion and insert the new one
                lc=open(inputFile,'w')
                for line in lines:
                    if 'MWEBV' in line:
                        lc.write('@MWEBV %5.4f\n' % (mwebv * extinctionFactor))
                    else:
                        lc.write(line)
                lc.close()
            else:
                print "WARNING: Zero or no extinction for %s" % (inputFile) 

            # Adjust the ZP reference, if required
            if not options.adjustZPref:
                continue

            adjust=None
            for line in lines:
                try:
                    entries=line.split()
                    if entries[0]=="@SURVEY":
                        if entries[1] in ['Riess1999_LC','Jha2006_LC']:
                            adjust=entries[1]
                            break
                except:
                    pass

            if adjust=='Riess1999_LC':
                cmd="sed 's/VEGA2/%s/' %s > temp1" % ('VEGA-R99',inputFile)
                
            elif adjust=='Jha2006_LC':
                cmd="sed 's/VEGA2/%s/' %s > temp1" % ('VEGA-J06',inputFile)

            if adjust!=None:
                sp.call(cmd,shell=True)
                cmd='cp temp1 %s' % (inputFile)
                sp.call(cmd,shell=True)

# We need to overwrite the file                

    return


if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="training.config",
                      help="Parameter file containting the location of various JLA parameters")

    parser.add_option("-a", "--adjustZPref", dest="adjustZPref", default=False,
                      action='store_true',
                      help="Adjust ZP reference")

    parser.add_option("-C", "--clean", dest="clean", default=False,
                      action='store_true',
                      help="Clean the directory")

    (options, args) = parser.parse_args()

    adjustExtinction(options)

