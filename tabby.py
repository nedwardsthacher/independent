# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:46:50 2017

@author: nickedwards
"""

# optimal aperture - line 880
# master bias - line 328
# secz - line 1044
# total flux - line 998

# ra dec -> xy: units degrees
# do_cal feed files, bias, flat, write out with tag _cal
# add comment to header for solved and cal and make a check file ()
# go back to optimal and fix it so it takes files

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

files, fct = tp.get_files(dir='/Users/nickedwards/python/independent/',prefix='KIC')

# 1. bias subtract
# subtract bias from each frame
def subtract_bias(files,dir='/Users/nickedwards/python/independent/'):
    # load bias
    # if goes into thacherphot need to add catch for no master_bias.fits
    bias_image, bias_header = fits.getdata(dir+'master_bias.fits',header=True)

    # go through all files and subtract bias from them
    data = {'image':np.zeros((len(files),2048,2048)), 'header':[]}
    for i in range(len(files)):
        image, header = fits.getdata(files[i],header=True)
        data['image'][i,:,:] = image - bias_image
        data['header'].append(header)

    # for divide_flat
    # maybe write new files after all calibration as _solved_calibrated.fits?
    return data

# 2. flat divide
# divide each image by the master flat
def divide_flat(files,dir='/Users/nickedwards/python/independent/',fil='V'):
    # load flat
    # if goes into thacherphot need to add catch for no master_flat.fits
    flat_image, flat_header = fits.getdata(dir+'master_flat_'+fil+'.fits', header= True)

    # take the files and bias subtract them
    #check to make sure the files and fil are the same?
    data = subtract_bias(files)

    for i in range(len(data['image'])): data['image'][i] /= flat_image

    return data
