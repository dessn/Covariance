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

import sys
import JLA_library as JLA

sys.path.append('../../../SpectroSN/MCMC/')
from cosmology import distance_modulus

def hubble_residuals(alpha, beta, M_1_B, sel):
    data = data_all[sel]
    zhel = data['zhel']
    mb = data['mb']
    x1 = data['x1']
    color = data['color']
    thirdvar = data['3rdvar']
    biascor = data['biascor']
    zcmb = data['zcmb']

    # observed distance modulus
    mu_obs = mb - biascor - (M_1_B - alpha * x1 + beta * color)

    # model distance modulus - assume Om_m = 0.3, w = -1??
    mu_theory = np.zeros_like(mu_obs)
    for j in range(len(mu_obs)):
        mu_theory[j] = distance_modulus(zhel[j], zcmb[j], 0.295, 0.705, -1, 0.0)
        # correct for host mass
        if thirdvar[j] > 10:
            mu_theory[j] = mu_theory[j] + Delta_M
    return mu_obs - mu_theory

def weights(C_eta, coh):
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

def reml(C_eta, coh_arr, sel):
    hr = hubble_residuals(alpha, beta, M_1_B, sel)
    reml = np.zeros_like(coh_arr)
    for j, coh in enumerate(coh_arr):
        w = weights(C_eta, coh)
        w=w[sel]
        ##w = w[w>0] ## temp: deal with negative weights??
        r = np.log(sum(w)) - sum([np.log(wi) for wi in w])
        for i, wi in enumerate(w):
            r += wi*hr[i]**2
        reml[j] = r
    return reml


def plot_reml(C_eta, sel, sel_label, ax):
    coh_arr = np.linspace(0.02, 0.20, 10)
    reml_arr = reml(C_eta, coh_arr, sel)
    ax.plot(coh_arr, reml_arr, 'm+', label=sel_label)
    ax.set_xlabel(r'$\sigma_i$')
    ax.set_ylabel('REML')
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
        else:
            sample.append('HST')

    return sample

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default=None,
                      help="Parameter file containing the location of various JLA parameters")

    (options, args) = parser.parse_args()

    # Read in the parameter file
    params = JLA.build_dictionary(options.config)

    # Use JLA values
    alpha = 0.14
    beta = 3.1
    M_1_B = -19.05
    Delta_M = -0.08
    # O_m?

    C_eta = fits.getdata(JLA.get_full_path(params['eta']))
    sigma_diag = np.genfromtxt(JLA.get_full_path(params['diag']), comments='#', \
        usecols=(0, 1, 2), dtype='f8,f8,f8', names=['coh', 'lens', 'pecvel'])

    # Read in the lightcurve fit parameters
    lcfits = JLA.get_full_path(params['lightCurveFits'])
    data_all = Table.read(lcfits)

    # Determine the survey from the target name
    data_all['sample']=get_sample(data_all['name'])

    zhel = data_all['zhel']
    names= data_all['name']
    sample=data_all['sample']

    print ' %d SNe total' % len(zhel)

    # split by sample

    sel_lowz1 = np.where((sample == 'low-z') & (zhel < 0.03))[0]
    sel_lowz2 = np.where((sample == 'low-z') & (zhel >= 0.03))[0]
    sel_SDSS1 = np.where((sample == 'SDSS') & (zhel < 0.2))[0]
    sel_SDSS2 = np.where((sample == 'SDSS') & (zhel >= 0.2))[0]
    sel_SNLS1 = np.where((sample == 'SNLS') & (zhel < 0.5))[0]
    sel_SNLS2 = np.where((sample == 'SNLS') & (zhel >= 0.5))[0]
    sel_HST  = np.where((sample == 'HST'))[0]
    selections = collections.OrderedDict()
    selections['low-z1']=sel_lowz1
    selections['low-z2']=sel_lowz2
    selections['SDSS1']=sel_SDSS1
    selections['SDSS2']=sel_SDSS2
    selections['SNLS1']=sel_SNLS1
    selections['SNLS2']=sel_SNLS2
    selections['HST']=sel_HST

    i = 0
    fig=plt.figure()
    gs = gridspec.GridSpec(4, 2)

    for sel_label, sel in selections.items():
        ax = fig.add_subplot(gs[i])
        plot_reml(C_eta, sel, sel_label, ax)
        i += 1
        ax.legend()

    plt.savefig('%s.png' % params['analysis'])
    plt.close()
