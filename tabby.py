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
from photutils.background import Background2D
from astropy import wcs
import pdb

def init():
    stars = np.array(['KIC8462852','SA41128','GD391E','GD391A','SA38326'])
    path = '../2017Mar16'
    files, fct = tp.get_files(dir=path,prefix='SA98-978')
    pdb.set_trace()
