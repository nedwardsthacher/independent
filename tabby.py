# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:46:50 2017

@author: nickedwards
"""

# optimal aperture - line 918


import matplotlib; matplotlib.use('Agg')
import thacherphot as tp
import thacherphot_ne as tpn
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropysics.coords import AngularCoordinate as angcor
from quick_image import *
import pdb
from datetime import datetime as dt
from matplotlib.dates import DateFormatter

stars = np.array(['KIC8462852','GD391E','GD391A','SA38326'])
date = np.array(['2017Mar23/','2017Apr02/','2017Apr05/','2017Apr06/'])

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

def reduc(files,dir='/Users/nickedwards/python/independent/'):

    # calibrate images
    cal = divide_flat(files,dir=dir)

    # huge dictionary for name, date, and the dict that optimal
    # aperture spits out, maybe append images and headers?
    data = {'name':[],'datetime':[],'optimal_aperture':[],
            'xcen':[],'ycen':[],'fwhm':[],'aspect':[],
            'snrmax':[],'totflux':[],'totflux_aperture':[],
            'chisq':[],'curve_of_growth':[],'secz':[]}

    # keys for optimal aperture dict but not any extra keys
    keys = np.array(data.keys())
    keys = keys[keys != 'datetime']
    keys = keys[keys != 'name']

    # go through the files and append name, date as date + UT, and
    # and optimal aperture indices to data
    for i in range(len(files)):
        ra = angcor(cal['header'][i]['RA']).d
        dec = angcor(cal['header'][i]['DEC']).d
        op = tpn.optimal_aperture(cal['image'][i],cal['header'][i],[20,30],ra=ra,dec=dec)
        for key in keys: data[key] = np.append(data[key],op[key])
        data['name'] = np.append(data['name'],cal['header'][i]['OBJECT'])
        datetime = cal['header'][i]['DATE'] + cal['header'][i]['UT']
        data['datetime'] = np.append(data['datetime'],datetime)

    return data

def plot_fluxes(data):

    plt.ion()
    plt.figure(1000)
    plt.clf()

    SA_i = data['name'] == 'SA38326'
    GDA_i = data['name'] == 'GD391A'
    GDE_i = data['name'] == 'GD391E'
    KIC_i = data['name'] == 'KIC8462852'

    SAflux = data['totflux'][SA_i]
    GDAflux = data['totflux'][GDA_i]
    GDEflux = data['totflux'][GDE_i]
    KICflux = data['totflux'][KIC_i]

    SAsecz = data['secz'][SA_i]
    GDAsecz = data['secz'][GDA_i]
    GDEsecz = data['secz'][GDE_i]
    KICsecz = data['secz'][KIC_i]

    SAline, = plt.plot(SAsecz,SAflux,'ro',label='SA 38-326')
    GDAline, = plt.plot(GDAsecz,GDAflux,'go',label='GD 391 A')
    GDEline, = plt.plot(GDEsecz,GDEflux,'yo',label='GD 391 E')
    KICline, = plt.plot(KICsecz,KICflux,'bo',label='KIC 8462852')

    plt.ylabel('totflux')
    plt.xlabel('secz')
    ax = plt.gca()
    ax.grid(True)
    plt.title('secz vs total flux on ' + data['datetime'][0][0:8])
    plt.legend(handles=[SAline,GDAline,GDEline,KICline],loc='best')

    plt.show()

    plt.figure(1001)
    plt.clf()

    SAtime = data['datetime'][SA_i]
    GDAtime = data['datetime'][GDA_i]
    GDEtime = data['datetime'][GDE_i]
    KICtime = data['datetime'][KIC_i]

    SAtime = [dt.strptime(time,'%d/%m/%y%H:%M:%S') for time in SAtime]
    GDAtime = [dt.strptime(time,'%d/%m/%y%H:%M:%S') for time in GDAtime]
    GDEtime = [dt.strptime(time,'%d/%m/%y%H:%M:%S') for time in GDEtime]
    KICtime = [dt.strptime(time,'%d/%m/%y%H:%M:%S') for time in KICtime]

    SAline, = plt.plot_date(SAtime,SAflux,'ro',label='SA 38-326')
    GDAline, = plt.plot_date(GDAtime,GDAflux,'go',label='GD 391 A')
    GDEline, = plt.plot_date(GDEtime,GDEflux,'yo',label='GD 391 E')
    KICline, = plt.plot_date(KICtime,KICflux,'bo',label='KIC 8462852')

    ax = plt.gca()
    xfmt = DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(xfmt)
    ax.grid(True)
    plt.ylabel('totflux')
    plt.xlabel('time (UT)')
    plt.title('secz vs total flux on ' + data['datetime'][0][0:8])
    plt.legend(handles=[SAline,GDAline,GDEline,KICline],loc='best')
    plt.show()

#files, fct = tp.get_files(dir='/Users/nickedwards/python/independent/',prefix='*',tag='solved')
#cal = divide_flat(files)
#data = reduc(files)
