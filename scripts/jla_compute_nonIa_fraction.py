"""Python program that computes the fraction of non-Ia in redshift bins
"""

from optparse import OptionParser

def getVariable(filename,variables):
    # Get the value of variables
    f=open(filename)
    lines=f.readlines()
    f.close()

    values=[]
    for line in lines:
        for variable in variables:
            entries=line.split()
            if entries[0]=='@'+variable:
                values.append(entries[1])

    return values

def compute_nonIa_fraction(options):

    import numpy
    import astropy.io.fits as fits
    from astropy.table import Table
    import JLA_library as JLA

    # -----------  Read in the configuration file ------------

    params = JLA.build_dictionary(options.config)

    # -----------  Read in the SNe ---------------------------

    SNeList = numpy.genfromtxt(options.SNlist,
                               usecols=(0, 2),
                               dtype='S30,S200',
                               names=['id', 'lc'])

    for i, SN in enumerate(SNeList):
        SNeList['id'][i] = SNeList['id'][i].replace('lc-', '').replace('.list', '')


    redshift=[]
    SNtype=[]
    for SN in SNeList:
        if 'DES' in SN['lc']:
            z,t=getVariable(SN['lc'],['Z_HELIO','SNTYPE'])
            redshift.append(float(z))
            SNtype.append(int(t))

    nSNe=len(redshift)
    types=numpy.zeros(nSNe,dtype=[('redshift','f4'),('SNtype','i4')])
    types['redshift']=redshift
    types['SNtype']=SNtype

    print "There are %d SNe" % (nSNe)



    bins=numpy.array([0.00,0.10,0.26,0.41,0.57,0.72,0.89,1.04])
    bias=numpy.array([0.00,0.015,0.024,0.024,0.024,0.023,0.026,0.025])
    selection=(types['SNtype']==1)
    SNeIa=numpy.histogram(types[selection]['redshift'],bins)
    SNeIa_total=numpy.histogram(types['redshift'],bins)

    if nSNe > SNeIa_total[0].sum():
        print 'Oops'
        exit()

    # Write out as an ascii table
    output=open(JLA.get_full_path(params['classification']),'w')
    output.write('# DES SN classification\n')
    output.write('# Written on %s\n' % (JLA.get_date()))
    output.write('# redshift\tRaw_bias\tfraction\n')
    output.write('%4.2f\t%4.3f\t%4.3f\n' % (0.0,0.0,0.0))
    for i in range(len(bins)-1):
        if SNeIa_total[0][i] > 0:
            output.write('%4.2f\t%4.3f\t%4.3f\n' % (bins[i+1],bias[i+1],1.0-1.0*SNeIa[0][i]/SNeIa_total[0][i]))

    output.close()
    

    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="JLA.config",
                      help="Parameter file containting the location of various JLA files")

    parser.add_option("-s", "--SNlist", dest="SNlist",
                      help="List of SN")


    (options, args) = parser.parse_args()


    compute_nonIa_fraction(options)
