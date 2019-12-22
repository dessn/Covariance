
"""Python program to compute the C_cal matrices
"""

from optparse import OptionParser
import subprocess as sp
import JLA_library as JLA
import numpy as np
import os
import shutil as sh
import time

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
    
    # Create the directory where the training will be done
    mkdir(params['trainingDir'])

    # Copy accross files from the initDir to the trainingDir
    # Consists of
    # i) Starting point for the mininimisation routine - pca_1_opt1_final.list
    # ii) Configuration files for pcafit - with and without error snake
    for file in os.listdir(params['initDir']):
        sh.copy(params['initDir']+'/'+file,params['trainingDir'])

    # Make the output directory
    # This is where we store the output surface.
    # It's name consists of the training sample used and the version of snpca used to do the training
    date=date=JLA.get_date()
    outputDir="/%s/data_%s_%s_%s/" % (params['outputDir'],date,params['trainingSample'],params['snpcaVersion'])
    mkdir(outputDir)

    # Change to the directory where the 
    os.chdir(params['trainingDir'])

    # Set the SALTPATH    
    os.environ['SALTPATH']=params['SALTPATH']
    print('SALT PATH is %s' %os.environ['SALTPATH'])
    
    # Part a) First training, withiout error snake

    # Step 1 - Train without the error snake
    #
    # Configuration parameter specified in training_without_error_snake.conf
    #
    # Start with an initial guess from pca_1_opt1_final.list, which has
    #   i) 1400 BSpline coefficents in M0: salt2_template_0.dat
    #      14 phases between -20 and +50, and 100 spline knots per phase = 14 x 100 =1400
    #  ii) 1400 BSpline coefficents in M1: salt2_template_1.dat
    #      14 phases between -20 and +50, and 100 spline knots per phase = 14 x 100 =1400
    # iii) 4 polynomial coefficients in colour law: salt2_color_correction.dat
    #      
    # Training list is provided in params['trainingList']

    cmd=['pcafit',
         '-l',params['trainingList'],
         '-c','training_without_error_snake.conf',
         '-p','pca_1_opt1_final.list',
         '-d']

    # Will produce a SALT2 surface, inlcuding
    # salt2_color_correction.dat - CL(lambda)
    # salt2_template_0.dat - M0
    # salt2_template_1.dat - M1
    # X0 * ( M0 + X1 * M1) * exp (C * CL(lambda))
    
    sp.call(' '.join(cmd),shell=True)

    
    # Step 2 - Compute uncertainties
    #
    # Uses the following output from step 1
    #
    # full_weight_1.fits is a 5956 x 5956 matrix made up of

    #   1400 (knots comp 0) + 2 (constraints on comp) + 1400 (knots comp 1) +
    #   2 (constraints on comp) + 4 (colour law) + 420 SNe (coord 0)  +
    #   420 SNe (coord 1) + 2 (constrains on coordinates) + 420 SNe (colours)
    #   + 1 886 (spectral recalibration)  = 5956
    
    #   index for comp. 0  [0,1399]
    #   index for comp. 1  [1402,2801]
    #   index for color law [2804,2807]
    #   index for Xi [2808,3647]
    #   index for color [3650,4069]
    #   index for spec. calib [4070,5955]

    # pca_1_opt1_final.list
    
    cmd=['write_pca_uncertainties', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits', 
         '2', 
         '1.0',
         '1.0']

    # This produces the following files
    # covmat_1_with_constraints.fits
    # color_law_covmat.dat
    # salt2_lc_relative_variance_0.dat
    # salt2_lc_relative_covariance_01.dat
    # salt2_lc_relative_variance_1.dat
    # salt2_spec_variance_0.dat
    # salt2_spec_covariance_01.dat
    # salt2_spec_variance_1.dat
    
    sp.call(' '.join(cmd),shell=True)

    # Step 3 - Compute error snake
    #
    # Uses the followng output from step 1
    #   full_weight_1.fits
    #   pca_1_opt1_final.list
    #
    # Uses the followng output from step 2
    #   covmat_1_with_constraints.fits is a 5956 x 5956 matrix
    
    cmd=['compute_error_snake',
         params['trainingList'],
         'training_without_error_snake.conf', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits',
         'covmat_1_with_constraints.fits']

    # model_covmat_for_error_snake.fits
    # spec_coverage.dat
    # spec_coverage_basis.fits
    # salt2_lc_dispersion_scaling.dat

    
    sp.call(' '.join(cmd),shell=True)

    sh.copy('pca_1_opt1_final.list', 'pca_1_opt1_final_first.list')
    sh.copy('model_covmat_for_error_snake.fits','model_covmat_for_error_snake_first.fits')
    sh.copy('salt2_lc_dispersion_scaling.dat', 'salt2_lc_dispersion_scaling_first.dat')
    
    # Part b Second training, with the error snake

    # The revised configration file contains two additional parameters
    # @ERRORSNAKE_COV_MAT model_covmat_for_error_snake_first.fits
    # @ERRORSNAKE_COV_SCALING salt2_lc_dispersion_scaling_first.dat
    
    # Step 4 - Second training using the output from the first three steps
    #
    # Uses the followng output from step 1
    #   pca_1_opt1_final_first.list
    #
    # Uses the followng output from step 3
    #    model_covmat_for_error_snake_first.fits 2800 x 2800 matix
    #    salt2_lc_dispersion_scaling_first.dat

    
    cmd=['pcafit',
         '-l',params['trainingList'],
         '-c','training_with_error_snake.conf',
         '-p','pca_1_opt1_final_first.list',
         '-d']
         
    sp.call(' '.join(cmd),shell=True)


    # Produces
    # salt2_template_0.dat
    # salt2_template_1.dat
    # salt2_color_correction.dat

    
    # Step 5 - Recompute uncertainties
    #
    # Uses the following output from step 4
    #
    # full_weight_1.fits
    # pca_1_opt1_final.list

    cmd=['write_pca_uncertainties', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits', 
         '2', 
         '1.0',
         '1.0']

    
    sp.call(' '.join(cmd),shell=True)

    # Step 6 - Recompute error snake
    #
    # Uses the followng output from step 4
    #   full_weight_1.fits
    #   pca_1_opt1_final.list
    #
    # Uses the followng output from step 5
    #   covmat_1_with_constraints.fits

    cmd=['compute_error_snake',
         params['trainingList'],
         'training_with_error_snake.conf', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits',
         'covmat_1_with_constraints.fits']

    sp.call(' '.join(cmd),shell=True)

    # Step 7 - Print out light curve residuals
    # 
    # Uses the followng output from step 4
    # full_weight_1.fits
    # pca_1_opt1_final.list
    # 
    # Uses the followng output from step
    # covmat_1_with_constraints.fits
    #
    
    cmd=['lightcurve_residuals',
         params['trainingList'],
         'training_with_error_snake.conf', 
         'pca_1_opt1_final.list',
         'full_weight_1.fits',
         'covmat_1_with_constraints.fits']
         
    sp.call(' '.join(cmd),shell=True)

    # Produces
    # lcresiduals.list

    
    # Step 8 - Post colour law fit
    #
    # Uses the following output from step 4
    # salt2_color_correction.dat
    # 
    # Uses the followng output from step 7
    # lcresiduals.list
    #
    # Produces seceral files, the most important of which is cle_final.list which is copied to salt2_color_dispersion.dat
    
    cmd=['post_color_law_fit',
         'lcresiduals.list', 
         'salt2_color_correction.dat',
         '-n',
         '3']

    sp.call(' '.join(cmd),shell=True)
    
    # Create the new surface
    for directory in ['Instruments','MagSys']:        
        sh.copytree(params['SALTPATH']+'/'+directory,outputDir+directory)

    for file in ['fitmodel.card']:
        sh.copy(params['SALTPATH']+'/'+file,outputDir+file)

    os.mkdir(outputDir+'salt2-4')

    
    # We seem to be missing the file salt2_color_dispersion.dat

    
    inputFiles=['salt2_color_correction.dat',
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

    for inputFile,outputFile in zip(inputFiles,outputFiles):
        sh.copy('%s' % (inputFile), "%s/salt2-4/%s" % (outputDir, outputFile))



    return

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-c", "--config", dest="config", default="SALT2.config",
                      help="configuration file containing SALT parameters")

    (options, args) = parser.parse_args()


    train_SALT2(options)
