######################################################
##Author:  Kellie McGuire     kellie@kelliejensen.com
##
##Takes FITS files and computes median value for each
##pixel and median value after eliminating outliers.
##Saves results to .txt files.
######################################################


from glob import glob
import numpy as np
from astropy.io import fits
import os.path
from os import path
import sys
import matplotlib.pyplot as plt
import time
import datetime
import argparse


###Create array of the data from FITS file
def get_data(fitsfiles):
    for i in range(len(fitsfiles)):
        hdulist = fits.open(fitsfiles[i])
        hdulist.info()
        fitsfiles[i] = hdulist[0].data
    fitsfiles = np.array(fitsfiles)
    hdulist.close()
    return fitsfiles #3d array

###Flattens 3d array into 2d
def flatten_files(fitsfiles):
    n_files, width, height = fitsfiles.shape
    files = np.zeros((n_files, width*height))
    for n in range(n_files):
        files[n] = fitsfiles[n].flatten()
    return files #2d array


#########################################################
###Calculates the median and stddev of each pixel over
###all files. Creates arrays populated with median
###and median after omitting outliers and saves to .txt
#########################################################
def median_imgs(flatfiles, width, height, outdir):
    n_files, file_length = flatfiles.shape


    # Vectorized version

    print("Make median image")
    medpix = np.median(flatfiles, axis=0)
    stdpix = np.std(flatfiles, axis=0)

    goodpixels = np.logical_and( flatfiles > (medpix - stdpix), flatfiles < (medpix + stdpix) )

    fullMedian = np.tile(medpix, (n_files, 1))
    print(np.sum(np.logical_not(goodpixels)))



    median_img = medpix

    print("Make median no outliers image")
    flatfilesWithoutOutliers = flatfiles
    flatfilesWithoutOutliers[np.logical_not(goodpixels)] = fullMedian[np.logical_not(goodpixels)]
    median_sans_outliers = np.median(flatfilesWithoutOutliers, axis=0)
    # for cell in range(file_length):


    #     # cells = flatfiles[:,cell:cell+1]   ###cells is a slice of the n_files pixel at index cell

    #     # median = np.median(cells)
    #     # stddev = np.std(cells)
    #     # filter_cells = []

    #     # ###remove outliers from cells
    #     # for element in cells:
    #     #     if element > median-stddev:
    #     #         if element < median+stddev:
    #     #             filter_cells.append(True)
    #     #         else: filter_cells.append(False)
    #     #     else:
    #     #         filter_cells.append(False)
    #     # cells_filtered = cells[filter_cells]

    #     # compute median of good cells only

    #     ###fill arrays
        
    #     median_sans_outliers[cell] = np.median(flatfiles[goodpixels[:, cell], cell])

    #     print("{} of {}".format(cell,file_length))

    median_img = median_img.reshape(width, height)
    median_sans_outliers = median_sans_outliers.reshape(width, height)
    fits.writeto(os.path.join(outdir, 'median_img.fits'), median_img)
    fits.writeto(os.path.join(outdir,'median_sans_outliers.fits'), median_sans_outliers,)
    return




def createMedianImages(filenames, outdir):


    start_time = time.time()


    data = get_data(filenames)

    #flatfiles = flatten_files(data)
    #n_files, width, height = data.shape

    #median_imgs(flatfiles,width,height, outdir)


    medianImage = np.median(data, axis=0)
    fits.writeto(os.path.join(outdir, "median-image.fits"), medianImage)


    print(f"time: {datetime.timedelta(seconds=int(time.time() - start_time))}")






if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Plot Skipper Image")
    parser.add_argument("-f", "--filenames", nargs="*", help="Skipper FITS file to plot")
    args = parser.parse_args()


    # distribute command line args to variables
    filenames = args.filenames

    outdir = os.path.split(filenames[0])[0]


    print(filenames)
    

    createMedianImages(filenames, outdir)
