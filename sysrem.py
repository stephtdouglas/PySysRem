# The following is an algorithm to remove systematic effects in a large set of
# light curves based on a paper by O. Tamuz, T. Mazeh and S. Zucker (2004)
# titled "Correcting systematic effects in a large set of photometric light
# curves".
#
# Originally written for PTF data by a student of Nick Law (we think), and
# modified by Marcel Agueros
#
# Now modified by Stephanie Douglas for K2 data

import scipy
import numpy as np
import sys
import math
import os.path
import os
import glob
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
	final_medians=np.zeros((stars_dim,2))

	nstars = stars_dim
	filenames = np.empty(nstars,"S100")
	filenames[:] = ""
	median_list = np.ones(nstars)*-9999

	# Import each of the star files
	for x,star in enumerate(star_list):

		# Remove flagged epochs from the light curve
		star.clean_up()

		final_medians[x,0] = star.median
		final_medians[x,1]=  star.std

		# Calculate residuals from the ORIGINAL light curve
		star_residuals = star.orig_mags - star.median

		# For the data points with quality flags != 0,
		# set the errors to a large value
		star_errors = star.orig_magerrors
		star_errors[star.good_mask==False] = 10**20

		# import the residual and error values into the matrices in the correct position (rows corresponding to stars, columns to epochs)
		residuals[x] = star_residuals
		errors[x] = star_errors

		# Save filenames into a list
		filenames[x] = star.filename.split("/")[-1]

	return residuals, errors, filenames

def sysrem(residuals, errors, filenames):

    stars_dim, epoch_dim = np.shape(residuals)

	# This medians.txt file is a 2D list with the first column being the medians of stars' magnitudes at different epochs (the good ones) and their standard deviations, so that they can be plotted against the results after errors are taken out below.
	num_errors=0
	while num_errors<5: 		# The number of linear systematics to remove
		c = np.zeros(stars_dim)
		a = np.zeros(epoch_dim)
		
		# minimize a and c values for a number of iterations, iter
		iter = 0
		while iter < 10:

			# Using the initial guesses for each a value of each epoch, minimize c for each star
			for star in range(stars_dim):
				numerator = 0
				denominator = 0
				for epoch in range(epoch_dim):
					numerator = numerator + (((a[epoch])*(residuals[star,epoch]))/((errors[star,epoch])**2))
					denominator = denominator + ((a[epoch])**(2))/((errors[star,epoch])**(2))
				c[star]=numerator/denominator

			# Using the c values found above, minimize a for each epoch
			for epoch in range(epoch_dim):
				numerator = 0
				denominator = 0
				for star in range(stars_dim):
					numerator = numerator + (((c[star])*(residuals[star,epoch]))/((errors[star,epoch])**2))
					denominator = denominator + ((c[star])**(2))/((errors[star,epoch])**(2))
				a[epoch]=numerator/denominator

			print iter, c[0], a[0]  # Write to the terminal one value of c and a to make sure things are running smoothly
			iter += 1


		# Create a matrix for the systematic errors:
		syserr=np.zeros((stars_dim, epoch_dim))
		for s in range(stars_dim):
			for e in range(epoch_dim):
				syserr[s,e]=a[e]*c[s]

		# Remove the systematic error
		residuals=residuals-syserr
		num_errors += 1

	outfile_performance = open('sysrem_performance.txt','w')

	# Reproduce the results in terms of medians and standard deviations for plotting
	perf=np.zeros((stars_dim, 2))
	x=0

	for filename in file_list:
		star = source_lc.source.from_ptf(filename)

		flag_list=star.flags
		good_residuals=[]
		final_median=[]
		std_dev=[]


		n = 0
		for flag in range(len(flag_list)):
			if flag_list[flag]==0:
				good_residuals.append(residuals[x][flag])
				star.mags[flag] = (good_residuals + median_list[x])[n]
				n+=1

		final_median = np.median(good_residuals + median_list[x])
		std_dev = np.std(good_residuals + median_list[x])
		print >> outfile_performance, final_median, std_dev
		perf[x,0]=final_median
		perf[x,1]=std_dev
		x += 1

		outf = open(filename.split(".")[0] + ".sysrem.txt","w")
		print >> outf, "# Source", filename.split("/")[-1].split("_")[1], "X:", star.xc, "Y:", star.yc
		mags = []
		for n in range(len(star.epochs)):
			print >> outf, "%.6f"%star.epochs[n],
			print >> outf, "%.5f"%star.mags[n],
			print >> outf, "%.5f"%star.magerrors[n],
			print >> outf, "%.7f"%star.xs[n],
			print >> outf, "%.7f"%star.ys[n],
			print >> outf, "%.7f"%star.xschip[n],
			print >> outf, "%.7f"%star.yschip[n],
			print >> outf, "%.7f"%star.flags[n]
		outf.close()

if __name__=="__main__":

	if len(sys.argv==0):
		print("Please provide a list of light curve files")
	elif len(sys.argv==1):
		# A filename has been provided, containing a list of files
		listfile = at.read(sys.argv[0])
		file_list = listfile["filename"]
	elif len(sys.argv>1):
		# Create a list of all the input files:
		file_list = sorted(sys.argv[1:])

	star_list = []
	for filename in file_list:
		star = source_lc.source.from_ptf(filename)
		star_list.append(star)

	residuals, errors, filenames = generate_matrix(star_list)
