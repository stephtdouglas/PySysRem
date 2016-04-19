# The following is an algorithm to remove systematic effects in a large set of
# light curves based on a paper by O. Tamuz, T. Mazeh and S. Zucker (2004)
# titled "Correcting systematic effects in a large set of photometric light
# curves".
#
# Originally written for PTF data by a student of Nick Law (we think), and
# modified by Marcel Agueros
#
# Now modified by Stephanie Douglas for K2 data

import sys
import os.path
import os

import astropy.io.ascii as at
import numpy as np

import source_lc

# Set the output print options:
np.set_printoptions(threshold=np.nan, precision=6)

def generate_matrix(star_list):

    stars_dim = len(star_list)
    epoch_dim = len(star_list[0].epochs)

    # Create empty matrices for residuals and corresponding errors with the found dimensions such that number of rows correspond to the number of available stars, and the number of columns correspond to each specific epoch:
    residuals = np.zeros((stars_dim, epoch_dim))
    errors = np.zeros((stars_dim, epoch_dim))
    medians = np.zeros(stars_dim)

    nstars = stars_dim
    filenames = np.empty(nstars,"S100")
    filenames[:] = ""
    median_list = np.ones(nstars)

    # Import each of the star files
    for x,star in enumerate(star_list):

        # Remove flagged epochs from the light curve
        star.clean_up()

        median_list[x] = star.median

        # Calculate residuals from the ORIGINAL light curve
        star_residuals = star.orig_mags - star.median

        # For the data points with quality flags != 0,
        # set the errors to a large value
        star_errors = np.copy(star.orig_magerrors)
        star_errors[star.good_mask==False] = 10**20

        # import the residual and error values into the matrices in the correct position (rows corresponding to stars, columns to epochs)
        residuals[x] = star_residuals
        errors[x] = star_errors

        # Save filenames into a list
        filenames[x] = star.filename.split("/")[-1]

    return residuals, errors, median_list, star_list

def sysrem(input_star_list):

    residuals, errors, median_list, star_list = generate_matrix(input_star_list)
    stars_dim, epoch_dim = np.shape(residuals)

    # This medians.txt file is a 2D list with the first column being the medians
    # of stars' magnitudes at different epochs (the good ones) and their
    # standard deviations, so that they can be plotted against the results after
    # errors are taken out below.
    for num_errors in range(5):         # The number of linear systematics to remove
        c = np.zeros(stars_dim)
        a = np.ones(epoch_dim)

        # minimize a and c values for a number of iterations, iter
        for iter in range(10):

            # Using the initial guesses for each a value of each epoch, minimize c for each star
            for s in range(stars_dim):
                err_squared = errors[s]**2
                numerator = np.sum(a*residuals[s]/err_squared)
                denominator = np.sum(a**2/err_squared)
                c[s] = numerator / denominator

            # Using the c values found above, minimize a for each epoch
            for ep in range(epoch_dim):
                err_squared = errors[:,ep]**2
                numerator = np.sum(c*residuals[:,ep]/err_squared)
                denominator = np.sum(c**2/err_squared)
                a[ep] = numerator / denominator

        # Create a matrix for the systematic errors:
        syserr=np.zeros((stars_dim, epoch_dim))
        for s in range(stars_dim):
            for e in range(epoch_dim):
                syserr[s,e]=a[e]*c[s]

        # Remove the systematic error
        residuals = residuals - syserr

    # Reproduce the results in terms of medians and standard deviations for plotting
    outfile_performance = open('sysrem_performance.txt','w')

    for x,star in enumerate(star_list):

        good_residuals = residuals[x][star.good_mask]
        correction = good_residuals + median_list[x]
        corrected_mags = np.copy(star.orig_mags)
        corrected_mags[star.good_mask] = correction

        final_median = np.median(good_residuals + median_list[x])
        std_dev = np.std(good_residuals + median_list[x])
        print >> outfile_performance, final_median, std_dev

        outfile_name = star.filename[:-4]+".sysrem.txt"
        data = {"epochs":star.orig_epochs,
                "mags":corrected_mags,
                "errors":star.orig_magerrors,
                "flags":star.orig_flags}
        formats = {"epochs":"%06f",
                "mags":"%0.5f",
                "errors":"%0.6f"}
        at.write(data,outfile_name,
                 names=["epochs","mags","errors","flags"],
                 formats=formats)

if __name__=="__main__":


    if len(sys.argv)==1:
        print("Please provide a list of light curve files")
    elif len(sys.argv)==2:
        # A filename has been provided, containing a list of files
        listfile = at.read(sys.argv[0])
        file_list = listfile["filename"]
    elif len(sys.argv)>2:
        # Create a list of all the input files:
        file_list = sorted(sys.argv[1:])

    star_list = []
    for filename in file_list:
        star = source_lc.source.from_ptf(filename)
        star_list.append(star)

    sysrem(star_list)
