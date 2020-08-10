"""Python wrapper to compare lightcurve parameters
"""

from optparse import OptionParser
import JLA_library as JLA
import numpy as np
import time
import os
import matplotlib.pyplot as plt

lcp=['DayMax','Rv','Tau','Color','X0','X1']
fitResults=["@CHI2_LC","@NDF_LC"]


class lightCurveFit:
    """A light curve fit"""
    def __init__(self,filename):
        f=open(filename)
        listing=f.readlines()
        f.close()

        self.lcp={}
        self.fitResults={}
        
        start=False
        for line in listing:
            entries=line.split()
            try:
                if entries[0]=="BEGIN_OF_FITPARAMS" and entries[1]=="Salt2Model":
                    start=True
                elif entries[0]=="END_OF_FITPARAMS" and entries[1]=="Salt2Model":
                    start=False
            except:
                pass
               
            if start:
                if entries[0] in lcp:
                    self.lcp[entries[0]]=float(entries[1])
                    self.lcp["%s_err" % (entries[0])]=float(entries[2])

            try:
                if entries[0] in fitResults:
                    self.fitResults[entries[0]]=float(entries[1])
            except:
                pass
                    
        return

def analyseLightCurveFits(options):

    # Get today's date
    date=date=JLA.get_date()

    # Read in the configuration file
    params=JLA.build_dictionary(options.config)

    # collate the lightcurves in the first directory and match with the light Curves in the second directory
    SNe={}

    dir1=os.listdir(params['dir1'])
    dir2=os.listdir(params['dir2'])

    for lcf in dir1:
        if lcf in dir2 and "sn" in lcf:
            SNname=lcf.replace(".dat","")
            SNe[SNname]=[lightCurveFit("%s/%s" % (params['dir1'],lcf)),lightCurveFit("%s/%s" % (params['dir2'],lcf))]
        else:
            print("Found %s in %s but not %s" % (lcf,params['dir1'],params['dir2']))

    # convert this into a function
    colour1=[]
    colour2=[]
    tau=[]
    Rv=[]
    for SN in SNe.keys():
        colour1.append(SNe[SN][0].lcp["Color"])
        colour2.append(SNe[SN][1].lcp["Color"])
        tau.append(SNe[SN][1].lcp["Tau"])
        Rv.append(SNe[SN][1].lcp["Rv"])


    # Plot the results
    fig=plt.figure(figsize=(12, 6))
    gs = fig.add_gridspec(2, 3)
    
    colmin=min(colour1+colour2)
    colmax=max(colour1+colour2)

    # Histogram above
    ax1=fig.add_subplot(gs[0,0])
    bins=np.arange(colmin,colmax+0.1,0.1)
    ax1.hist(colour2,bins,alpha=0.5,label="With",color='b')
    ax1.hist(colour1,bins,alpha=0.5,label="Without",color='orange')
    y_min, y_max = ax1.get_ylim()
    ax1.plot([np.median(colour2),np.median(colour2)],[0,y_max],color='b',ls='dashed')
    ax1.plot([np.median(colour1),np.median(colour1)],[0,y_max],color='orange',ls='dashed')
    ax1.legend()

    # Plot below
    ax2=fig.add_subplot(gs[1,0])
    ax2.plot(colour1,colour2,'b.')
    ax2.plot([colmin,colmax],[colmin,colmax])
    ax2.set_xlabel("Colour - without host galaxy extinction")
    ax2.set_ylabel("Colour - with\n host galaxy extinction")

    # Histogram above
    ax3=fig.add_subplot(gs[0,1])
    taumin=min(tau)
    taumax=max(tau)
    bins=np.arange(taumin,taumax+0.1,0.1)
    ax3.hist(tau,bins,alpha=0.5,label="With",color='b')
    y_min, y_max = ax3.get_ylim()
    ax3.fill_between([taumin,0],[y_max,y_max],hatch='/',alpha=0.5,color='r')
    ax3.legend()

    # Plot below
    ax4=fig.add_subplot(gs[1,1])
    ax4.plot(tau,colour2,'b.')
    ax4.set_xlabel("E(B-V) - with host galaxy extinction")

    # Histogram above
    ax5=fig.add_subplot(gs[0,2])
    Rvmin=min(Rv)
    Rvmax=max(Rv)
    bins=np.arange(Rvmin,Rvmax+0.1,0.1)
    ax5.hist(Rv,bins,alpha=0.5,label="With",color='b')
    y_min, y_max = ax5.get_ylim()
    ax5.plot([np.median(Rv),np.median(Rv)],[0,y_max],color='b',ls='dashed')
    ax5.plot([3.1,3.1],[0,y_max],color='orange',ls='dashed')
    ax5.set_xlabel("R_V - with\n host galaxy extinction")
    
    
    ax5.legend()
    
    plt.show()
    plt.close()

    return



if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="SNFIT.config",
                      help="configuration file containing SALT parameters")

    (options, args) = parser.parse_args()


    analyseLightCurveFits(options)

