import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import argparse
import sys
import palettable

sys.path.append("/home/apiers/damicm/damicm-image-preproc/AutoAnalysis")

cmap = palettable.cmocean.sequential.Thermal_20.mpl_colormap
import readFits

halfPixArr = 3081

if __name__ == '__main__':

	

	parser = argparse.ArgumentParser(description="Plot Skipper Image")
	parser.add_argument("-f", "--filename", help="Skipper FITS file to plot")
	args = parser.parse_args()

	# distribute command line args to variables
	filename = args.filename

	header, data = readFits.read(filename)

	data = data[:, :, 0]
	nrows, ncols = data.shape

	fig, ax = plt.subplots(1, 1, figsize=(12, 8))

	ax.set_xlabel("x (pixels)", fontsize=18)
	ax.set_ylabel("y (pixels)", fontsize=18)
	ax.tick_params(labelsize=16)

	# Perform median subtraction of both overscans
	data[:, :ncols//2] -= np.median(data[:, (ncols//2 - 100):ncols//2])
	data[:, ncols//2:] -= np.median(data[:, (ncols//2):(ncols//2 + 100)])


	# physicalData = np.zeros(data.shape)
	# physicalData[:, :ncols//2] = np.fliplr(data[:, :ncols//2]) # L amplifier
	# physicalData[:, ncols//2:] = np.fliplr(data[:, ncols//2:]) # U amplifier

	# Plot median image
	cax = ax.imshow(data[5:,:], vmin=-60, vmax=20, origin="lower", cmap=cmap)
	cb = fig.colorbar(cax, ax=ax)
	cb.set_label("Pixel Value (ADU)", fontsize=16)


	plt.show()
