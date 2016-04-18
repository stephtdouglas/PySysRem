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
import numpy
import sys
import math
import os.path
import os
import glob
import source_lc


# Set the output print options:
numpy.set_printoptions(threshold=numpy.nan, precision=6)

# Create a list of all the input files:
file_list = sorted(sys.argv[1:])

one_source = source_lc.source(file_list[0]) 		# import one of the files
epoch_list1 = one_source.epochs 			# make a list of the epochs in that file
epoch_dim = len(epoch_list1)				# length of the list of epochs

stars_dim = len(file_list)	# find the number of files available in a given directory(quad) which corresponds to the number of stars available

# Create empty matrices for residuals and corresponding errors with the found dimensions such that number of rows correspond to the number of available stars, and the number of columns correspond to each specific epoch:
residuals = numpy.zeros((stars_dim, epoch_dim))
errors = numpy.zeros((stars_dim, epoch_dim))

# Initialize specific parameters:
x = 0
final_medians=numpy.zeros((stars_dim,2))
stars=[]
median_list=[]

# Import each of the star files
for filename in file_list:
	star = source_lc.source(filename)		# Get the source
	mag_list=star.mags				# Get the magnitudes
	flag_list=star.flags				# Get the Quality Flags
	error_list = star.magerrors			# Get the errors
	good_magnitudes=[]				# Create an empty list for good magnitudes to be stored in
	std_dev=[]

	# If the Quality Flag is = 0 then import the corresponding magnitude in the good_magnitude list
	for flag in range(len(flag_list)):
		if flag_list[flag]==0:
			good_magnitudes.append(mag_list[flag])

	median = numpy.median(good_magnitudes)		# Find the median of good magnitudes
	median_list.append(median)			# Import each star's median in a list
	std_dev=numpy.std(good_magnitudes)

	final_medians[x,0]=median
	final_medians[x,1]=std_dev
	r = mag_list - median				# Residual = magnitude - median

	# For the data points with quality flags != 0, set the errors to a large value
	for item in range(len(flag_list)):
		if flag_list[item] != 0:
			error_list[item] = 10**20

	# import the residual and error values into the matrices in the correct position (rows corresponding to stars, columns to epochs)
	for y in range(epoch_dim):
		residuals[x,y] = r[y]
		errors[x,y] = error_list[y]
	x += 1

	# Import filenames into a list named stars, to keep track of what row corresponds to which star: the first element = first row of matrix and so on
	stars.append(filename)

# This medians.txt file is a 2D list with the first column being the medians of stars' magnitudes at different epochs (the good ones) and their standard deviations, so that they can be plotted against the results after errors are taken out below.
num_errors=0
while num_errors<5: 		# The number of linear systematics to remove
	c=[]
	a=[]

	# Create a list for star errors, with initial values of zero
	for i in range(stars_dim):
		c.append(0)

	# Create a list for epoch errors, with initial guess of 1 for all epochs
	for i in range(epoch_dim):
		a.append(1)

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
	syserr=numpy.zeros((stars_dim, epoch_dim))
	for s in range(stars_dim):
		for e in range(epoch_dim):
			syserr[s,e]=a[e]*c[s]

	# Remove the systematic error
	residuals=residuals-syserr
	num_errors += 1

outfile_performance = open('sysrem_performance.txt','w')

# Reproduce the results in terms of medians and standard deviations for plotting
perf=numpy.zeros((stars_dim, 2))
x=0

for filename in file_list:
	star = source_lc.source(filename)

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

	final_median = numpy.median(good_residuals + median_list[x])
	std_dev = numpy.std(good_residuals + median_list[x])
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
		print >> outf, "%.7f"%star.flags[n],
		print >> outf, star.images[n]
	outf.close()
