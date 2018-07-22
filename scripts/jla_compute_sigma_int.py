"""
"""

# from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.io.fits as fits
from astropy.table import Table
from optparse import OptionParser
import collections
import scipy.optimize as sp

import sys
import JLA_library as JLA

sys.path.append('../../../SpectroSN/MCMC/')
from cosmology import distance_modulus

def hubble_residuals(alpha, beta, M_1_B, sel, sel_label, ax2, C_eta, sigma_diag):
    data = data_all[sel]
    zhel = data['zhel']
    mb = data['mb']
    x1 = data['x1']
    color = data['color']
    thirdvar = data['3rdvar']
    biascor = data['biascor']
    zcmb = data['zcmb']

    N = len(C_eta)/3
    C_mu = np.zeros((N,N))
    for i in range(N):
        C_mu[i, i] += sigma_diag['lens'][i]**2 \
                   +  sigma_diag['pecvel'][i]**2
    A = JLA.A_fn(alpha,beta,N)
    C_mu += A.dot(C_eta.dot(A.T))
    error=np.sqrt(np.diag(C_mu))[sel]


    # observed distance modulus
    mu_obs = mb - biascor - (M_1_B - alpha * x1 + beta * color)

    # model distance modulus - assume Om_m = 0.3, w = -1??
    mu_theory = np.zeros_like(mu_obs)
    for j in range(len(mu_obs)):
        mu_theory[j] = distance_modulus(zhel[j], zcmb[j], 0.295, 0.705, -1, 0.0)
        # correct for host mass
        if thirdvar[j] > 10:
            mu_theory[j] = mu_theory[j] + Delta_M

    offset=np.average(mu_obs - mu_theory)
    ax2.errorbar(zhel, mu_obs - mu_theory - offset, error, marker='o', mfc='red', mec='red', ms=3, label=sel_label,ls='None')
    ax2.set_xlabel('redshift')
    ax2.set_ylim(-0.8,0.8)

    return mu_obs - mu_theory - offset

def weights(C_eta, coh, sigma_diag):
    N = len(C_eta)/3
    C_mu = np.zeros((N,N))
    for i in range(N):
        C_mu[i, i] += coh**2 + sigma_diag['lens'][i]**2 \
                     + sigma_diag['pecvel'][i]**2
    A = JLA.A_fn(alpha,beta,N)
    C_mu += A.dot(C_eta.dot(A.T))
    w = np.zeros(N)
    for i in range(N):
        w[i] = 1/C_mu[i][i]
    return w

def reml(C_eta, coh_arr, sel, sel_label, ax2, sigma_diag):
    hr = hubble_residuals(alpha, beta, M_1_B, sel, sel_label, ax2, C_eta, sigma_diag)
    reml = np.zeros_like(coh_arr)
    for j, coh in enumerate(coh_arr):
        w = weights(C_eta, coh, sigma_diag)
        w=w[sel]
        ##w = w[w>0] ## temp: deal with negative weights??
        r = np.log(sum(w)) - sum([np.log(wi) for wi in w])
        for i, wi in enumerate(w):
            r += wi*hr[i]**2
        reml[j] = r
    return reml


def plot_reml(C_eta, sel, sel_label, ax, ax2, sigma_diag):
    coh_arr = np.linspace(0.00, 0.30, 16)
    reml_arr = reml(C_eta, coh_arr, sel, sel_label, ax2, sigma_diag)
    ax.plot(coh_arr, reml_arr, 'm+', label=sel_label)
    ax.set_xlabel(r'$\sigma_i$')
    ax.set_ylabel('REML')
    ax.set_xlim(-0.02,0.30)

    # Fit a polynomial to the points
    coeff=np.polyfit(coh_arr, reml_arr, 4)
    p=np.poly1d(coeff)
    coh_arr2=np.linspace(0.00, 0.20, 201)
    fit=p(coh_arr2)
    ax.plot(coh_arr2,fit)

    # Determine where the minumum occurs
    minimum=sp.fmin(p,0.1,disp=False)
    print "Miminum at %4.3f" % (minimum[0])
    # The uncetainty occurs whwn the the REML increases by 1 - need to double check
    # To find where the the REML increases by 1, we look for the roots
    coeff=np.polyfit(coh_arr, reml_arr-p(minimum[0])-1, 4)
    p=np.poly1d(coeff)
    sol=sp.root(p,[0.0,0.2])

    m=minimum[0]
    upper=sol.x[1]-m
    lower=m-sol.x[0]
    ax.plot(coh_arr2,fit,label="Min at %4.2f+%4.2f-%4.2f" % (m,upper,lower))
    
    return

def get_sample(names):
    sample=[]
    for name in names:
        if 'SDSS' in name:
            sample.append('SDSS')
        elif 'sn' in name[0:2]:
            sample.append('low-z')
        elif name[2:4] in ['D1','D2','D3','D4']:
            sample.append('SNLS')
        elif 'DES' in name:
            sample.append('DES')
        else:
            sample.append('HST')

    return sample

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default=None,
                      help="Parameter file containing the location of various files")

    (options, args) = parser.parse_args()

    # Read in the parameter file
    params = JLA.build_dictionary(options.config)

    # Use JLA values
    alpha = 0.14
    beta = 3.1
    M_1_B = -19.02 # We use a brigter value JLA is -19.05
    Delta_M = -0.08
    # O_m?

    C_eta = fits.getdata(JLA.get_full_path(params['eta']))
    sigma_diag = np.genfromtxt(JLA.get_full_path(params['diag']), comments='#', \
        usecols=(0, 1, 2), dtype='f8,f8,f8', names=['coh', 'lens', 'pecvel'])

    # Read in the lightcurve fit parameters
    lcfits = JLA.get_full_path(params['lightCurveFits'])
    data_all = Table.read(lcfits)

    # Read in the list of SNe that are used in the BBC analysis
    include=np.genfromtxt(JLA.get_full_path(params['include']), comments='#', \
        usecols=(0), dtype=[('name','a20')])


    include2=[]
    for name in include['name']:
        if name[0]=='1':
            include2.append("DES_0%s" % name)
        else:
            include2.append("sn%s" % name)

    print "There are %d SNe in the DES3YR analysis" % (len(include2))

    # Read in the order in which the SNe appear in the files
    ordering=np.genfromtxt(JLA.get_full_path(params['order']), comments='#', \
        usecols=(0), dtype=[('name','a30')])

    
    nSNe=len(ordering)
    ##tempList=[]
    
    keep = np.zeros(3*nSNe,bool)
    # Modify the input, so that only the SNe that pass cuts are used ...
    for i,SN in enumerate(ordering['name']):
        if SN.replace('lc-','').replace('.list','') in include2:
            keep[3*i:3*i+3]=True
            ##tempList.append(SN.replace('lc-','').replace('.list',''))

    # The light curve fits file has 561 SNe, so we do that one separately
    # This will need to be updated by Tamara. We can then remove this code
    keepfits=np.zeros(len(data_all),bool)
    # Modify the input, so that only the SNe that pass cuts are used ...
    for i,SN in enumerate(data_all['name']):
        if SN in include2:## and SN in tempList:
            keepfits[i]=True


    print keep.sum(),keepfits.sum()

    data_all=data_all[keepfits]

    # Determine the survey from the target name
    data_all['sample']=get_sample(data_all['name'])

    zhel = data_all['zhel']
    names= data_all['name']
    sample=data_all['sample']

    print ' %d SNe total' % len(zhel)

    # split by sample

    sel_lowz1 = np.where((sample == 'low-z') & (zhel < 0.03))[0]
    sel_lowz2 = np.where((sample == 'low-z') & (zhel >= 0.03))[0]
#    sel_SDSS1 = np.where((sample == 'SDSS') & (zhel < 0.2))[0]
#    sel_SDSS2 = np.where((sample == 'SDSS') & (zhel >= 0.2))[0]
#    sel_SNLS1 = np.where((sample == 'SNLS') & (zhel < 0.5))[0]
#    sel_SNLS2 = np.where((sample == 'SNLS') & (zhel >= 0.5))[0]
    sel_DES1 = np.where((sample == 'DES') & (zhel < 0.5))[0]
    sel_DES2 = np.where((sample == 'DES') & (zhel >=0.5))[0]
##    sel_DES3 = np.where((sample == 'DES') & (zhel >= 0.5))[0]
#    sel_HST  = np.where((sample == 'HST'))[0]
    selections = collections.OrderedDict()
    selections['low-z1']=sel_lowz1
    selections['low-z2']=sel_lowz2
#    selections['SDSS1']=sel_SDSS1
#    selections['SDSS2']=sel_SDSS2
#    selections['SNLS1']=sel_SNLS1
#    selections['SNLS2']=sel_SNLS2
    selections['DES1']=sel_DES1
    selections['DES2']=sel_DES2
##    selections['DES3']=sel_DES3
#    selections['HST']=sel_HST

    i = 0
    fig=plt.figure()
    fig2=plt.figure()
    gs = gridspec.GridSpec(2, 2)
    gs2 = gridspec.GridSpec(2, 2)

    print len(sel_lowz1), len(sel_lowz2), len(sel_DES1), len(sel_DES2)##, len(sel_DES3)

    keepEta=np.matrix(keep).T*np.matrix(keep)
    nSNe=keep.sum()
    C_eta_sub=C_eta[keepEta].reshape(nSNe,nSNe)
    print C_eta_sub.shape

    for sel_label, sel in selections.items():
        ax = fig.add_subplot(gs[i])
        ax2 = fig2.add_subplot(gs2[i])
        plot_reml(C_eta_sub, sel, sel_label, ax, ax2, sigma_diag)
        i += 1
        ax.legend()

    fig.savefig('%s.png' % params['analysis'])
    fig2.savefig('%s_resid.png' % params['analysis'])
    plt.close()
