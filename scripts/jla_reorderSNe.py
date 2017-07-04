"""
Python program to reorder the SNe that appear in Cstat
"""

from optparse import OptionParser

def reorderSNe(options):
    # The ordering of the SNe produced by salt2_stat does not reflect the order they were input
    # The ordering in the output file is written to the file sne_mu.list
    # SDSS SNe in the JLA sample get called @SN 12856.0, which means that the output of salt2_stat has these names.
    # Some nearby SNe have sn in front of their names. Others do not 
    # In one case a lower case v is used
    # DES SNe in the DES sample are listed as 01248677 

    import numpy
    import astropy.io.fits as fits
    import JLA_library as JLA
    from astropy.table import Table

    # -----------  Read in the SN ordering ------------------------
    SNeList = Table(numpy.genfromtxt(options.SNlist,
                               usecols=(0, 2),
                               dtype='S30,S200',
                               names=['id', 'lc']))


    # -----------  Read in the file that species the ordering of the matrix produced by Cstat ------------
    statList = Table(numpy.genfromtxt(options.input,
                                      usecols=(0,1),
                                      dtype='S30,float',
                                      names=['id','z'],skip_header=13))


    # We use the -9 as a way to catch errors.
    reindex=numpy.zeros(len(statList),int)-9

    for i,SNname in enumerate(SNeList['id']):
        name=SNname.replace("lc-","").replace(".list","")
        for j,SNname2 in enumerate(statList['id']):
            if SNname2 == name or ("SDSS"+SNname2.replace(".0","") == name and "SDSS" in name) or \
                    (SNname2.replace("sn","")==name.replace("sn","")) or ("DES" in name and "DES_0"+SNname2==name):
                if reindex[i]!=-9:
                    print SNname,SNname2
                reindex[i]=j
                #print SNname,SNname2,i,j
                
    for index,value in enumerate(reindex):
        if value==-9:
            print "Error"
            print index,SNeList[index]['id']
            exit()

    print "The numbers should be the same"
    print len(SNeList), len(statList), len(numpy.unique(reindex))

    # -----------  Read in Cstat and re-order --------------------------------------

    Cstat=fits.getdata(options.file)

    # We use brute force to reorder the elements
    # Recall that for each SNe, there is an error associated with the peak mag, colour and stretch

    Cstat_new=numpy.copy(Cstat) * 0.0

    nSNe=len(SNeList)
    for i in  range(nSNe):
        for j in range(nSNe):
            Cstat_new[3*reindex[i]:3*reindex[i]+3,3*reindex[j]:3*reindex[j]+3]=Cstat[3*i:3*i+3,3*j:3*j+3]


    date = JLA.get_date()
    fits.writeto('%s_Cstat_%s.fits' % (options.prefix,date),Cstat_new,clobber=True) 

    return None

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-s", "--SNlist", dest="SNlist", default=None,
                      help="List of SN")

    parser.add_option("-p", "--prefix", dest="prefix", default="JLA++host",
                      help="Prefix to SN name")

    parser.add_option("-f", "--file", dest="file", default=None,
                      help="The input FITS file")

    parser.add_option("-i", "--input", dest="input", default=None,
                      help="The ordering as produced by salt2_stat")
    
    (options, args) = parser.parse_args()

    reorderSNe(options)
