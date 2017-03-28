# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:46:50 2017

@author: nickedwards
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import robust as rb
import thacherphot as tp
import hcongrid as h
from astropy import wcs
import pdb
import matplotlib; matplotlib.use('Agg')

stars = np.array(['KIC8462852','SA41128','GD391E','GD391A','SA38326'])
date = np.array(['../2017Mar16','../2017Mar23'])

# based on the dates and stars above
# astrometry() does the astrometry for all the stars for all the dates
#can import thacherphot as tp and just do it on bellerophon
# get files and do_astrometry(files)

def atrometry():
    for i in range(len(date)):
        for j in range(len(stars)):
            files, fct = tp.get_files(dir=date[i],prefix=stars[j])
            tp.do_astrometry(files,clobber=True)
        biasfiles, bfct = tp.get_files(dir=date[i],prefix='bias')

# send find flux the astrometrically solved images and
# it runs the tot_flux function from thacherphot
def find_flux(files):
    # slice file names into usful name then make a dict
    for file in range(len(files)):
        # appenf the dict based on files name
