# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:46:50 2017

@author: nickedwards
"""

# optimal aperture - line 880
# master bias - line 328
# secz - line 1044

import matplotlib; matplotlib.use('Agg')
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import robust as rb
import thacherphot as tp
import hcongrid as h
from astropy import wcs
import pdb

stars = np.array(['KIC8462852','SA41128','GD391E','GD391A','SA38326'])
date = np.array(['../2017Mar23'])

# 0. astrometrically solve images
# based on the dates and stars above
# can import thacherphot as tp and just do it on bellerophon
# get files and do_astrometry(files)

# 1. bias subtract
# subtract bias from each frame
def subtract_bias(files,dir='/Users/nickedwards/python/independent/'):
    # load bias
    # if goes into thacherphot need to add catch for no master_bias.fits
    bias_image, bias_header = fits.getdata(dir+'master_bias.fits',header=True)

    # go through all files and subtract bias from them
    # is there a better way of doing this?
    images = np.array([])
    headers = []
    for i in range(len(files)):
        image, header = fits.getdata(files[i],header=True)
        pdb.set_trace()
        image = image - bias_image
        # need to work on adding 2D array into this so that images[i] is an image
        images = np.append(images, image)
        headers = np.append(headers, header)

    # for divide_flat
    # maybe write new files after all calibration as _solved_calibrated.fits?
    return images, headers

# 2. flat divide
# divide each image by the master flat
def divide_flat(files,dir='/Users/nickedwards/python/independent/',fil='V'):
    # load flat
    # if goes into thacherphot need to add catch for no master_flat.fits
    flat_image, flat_header = fits.getdata(dir+'master_flat_'+fil+'.fits', header= True)

    # take the files and bias subtract them
    #check to make sure the files and fil are the same?
    images, headers = subtract_bias(files)

    for i in range(len(images)): images[i] /= flat_image

    return images, headers
