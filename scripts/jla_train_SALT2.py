"""Python program to compute the C_cal matrices
"""

from optparse import OptionParser
import subprocess as sp
import JLA_library as JLA
import numpy as np
import os
import shutil as sh

# Usage
# jla_train_SALT2.py optioons
# 
# The method follows the procedures specified in e-mail sent by M. Betoule on 15/11/2017

def mkdir(dir):
    try:
        os.mkdir(dir)
    except:
        pass

    return

def train_SALT2(options):
    # Read in the configuration file
    params=JLA.build_dictionary(options.config)
    
    # Make the initialisation and training directories
    mkdir(params['initDir'])
    mkdir(params['trainingDir'])

    # Copy accross files from the initDir to the trainingDir
    for file in os.listdir(params['initDir']):
        sh.copy(params['initDir']+'/'+file,params['trainingDir'])
    
    os.chdir(params['trainingDir'])

    # Part a) First training

    # Step 1 - Train without the error snake
    cmd=['pcafit',
         '-l','trainingsample_snls_sdss_v5.list',
         '-c','training_without_error_snake.conf',
         '-p','pca_1_opt1_final.list',
         '-d']

##    sp.call(' '.join(cmd),shell=True)

    # Step 2 - Compute uncertainties
    cmd=['write_pca_uncertainties', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits', 
         '2', 
         '1.0',
         '1.0']

##    sp.call(' '.join(cmd),shell=True)

    # Step 3 - Compute error snake
    cmd=['Compute_error_snake',
         'trainingsample_snls_sdss_v5.list', 
         'training_without_error_snake.conf', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits',
         'covmat_1_with_constraints.fits']
##    sp.call(' '.join(cmd),shell=True)

##    sh.copy('pca_1_opt1_final.list', 'pca_1_opt1_final_first.list')
##    sh.copy('model_covmat_for_error_snake.fits','model_covmat_for_error_snake_first.fits')
##    sh.copy('salt2_lc_dispersion_scaling.dat', 'salt2_lc_dispersion_scaling_first.dat')
    # Part b Second training, with the error snake

    # Step 4 - Second training using the output from the first three steps
    cmd=['pcafit',
         '-l','trainingsample_snls_sdss_v5.list',
         '-c','training_with_error_snake.conf',
         '-p','pca_1_opt1_final_first.list',
         '-d']
         
##    sp.call(' '.join(cmd),shell=True)

    # Step 5 - Recompute uncertainties
    cmd=['write_pca_uncertainties', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits', 
         '2', 
         '1.0',
         '1.0']

##    sp.call(' '.join(cmd),shell=True)

    # Step 6 - Recompute error snake
    cmd=['Compute_error_snake',
         'trainingsample_snls_sdss_v5.list', 
         'training_with_error_snake.conf', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits',
         'covmat_1_with_constraints.fits']

##    sp.call(' '.join(cmd),shell=True)

    # Step 7 - Print out light curve residuals
    cmd=['lightcurve_residuals',
         'trainingsample_snls_sdss_v5.list', 
         'training_with_error_snake.conf', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits',
         'covmat_1_with_constraints.fits']

##    sp.call(' '.join(cmd),shell=True)

    # Step 8 - Post colour law fit
    cmd=['post_color_law_fit',
         'lcresiduals.list', 
         'salt2_color_correction.dat',
         '-n',
         '-3']

##    sp.call(' '.join(cmd),shell=True)

    # Step 9 - Draw templates
    # We have our own code for this step.

    # Watch out for the name changes

    inputFiles=['salt2_color_correction_final.dat',
                'salt2_lc_relative_variance_0.dat',
                'salt2_spec_variance_1.dat',
                'cle_final.list',
                'salt2_lc_relative_variance_1.dat',
                'salt2_template_0.dat',
                'salt2_lc_dispersion_scaling.dat',
                'salt2_spec_covariance_01.dat',
                'salt2_template_1.dat',
                'salt2_lc_relative_covariance_01.dat',
                'salt2_spec_variance_0.dat']

    outputFiles=['salt2_color_correction.dat',
                'salt2_lc_relative_variance_0.dat',
                'salt2_spec_variance_1.dat',
                'salt2_color_dispersion.dat',
                'salt2_lc_relative_variance_1.dat',
                'salt2_template_0.dat',
                'salt2_lc_dispersion_scaling.dat',
                'salt2_spec_covariance_01.dat',
                'salt2_template_1.dat',
                'salt2_lc_relative_covariance_01.dat',
                'salt2_spec_variance_0.dat']

# We seem to be missing a file

    for inputFile,outputFile in zip(inputFiles,outputFiles):
        sh.copy('%s' % (inputFile), "%s/%s" % ('../snfit_data/salt2-4/', outputFile))



    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="SALT2.config",
                      help="configuration file containing SALT parameters")

    (options, args) = parser.parse_args()


    train_SALT2(options)
