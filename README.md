# Covariance

Code to compute the individual SN covariance matrices in JLA (Betoule+ 2014 Section 5.5) 

---

To run the scripts you need a configuration file and list of SNe in the Covariance/scripts/ directory. Here are a   few lines from each of mine (for the JLA sample) to indicate the format:

JLA.config:

# Config file for JLA sample
lightCurves   	      $JLA/data/jla_light_curves/     # the JLA lightcurves
lightCurveFits        $JLA/data/JLA_20160531.fits     # fits file containing results of lightcurve fits

where I have a $JLA environmental variable pointing to where these data are. For individual scripts you may need other keys in this configuration file - contact Bonnie or Chris if you are trying to run a script and something is missing.

JLA.list:

lc-03D1au.list LC $JLA/data/jla_light_curves_dm/new_lc-03D1au.list
lc-03D1aw.list LC $JLA/data/jla_light_curves_dm/new_lc-03D1aw.list
lc-03D1ax.list LC $JLA/data/jla_light_curves_dm/new_lc-03D1ax.list

---

Scripts:

- jla_compute_C*.py

code for C_* for each contribution e.g. host, dust, bias, nonia, model, cal, pecvel, stat

general usage is 'python jla_compute_C*.py -s JLA.list' (or other list of SNe). Some of these require other arguments and input. We're working on making these available; in the meantime email us. See docstring in each code for technical details of each calculation.

- jla_compute_diag_terms.py

computes diagonal terms in magnitude covariance matrix (note: all other covariance matrices are in eta i.e. m_B, X1, C) due to lensing, scatter in peculiar velocity and intrinsic scatter

- jla_compute_rel_size.py

computes relative contributions of various covariance matrices (each systematic and statistics) according to method in B14 Section 6.2

- JLA_library.py

library of common functions used in above scripts


The following scripts aren't essential parts of the chain:

- jla_merge_lightcurve_fits.py

formats lightcurves from multiple files and/or samples into .fits files

- jla_compute_DateofMax.py

inserts date of maximum into light curve files

- jla_compute_ZP.py

computes magnitude of CALSPEC standards