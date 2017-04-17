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

    keys = [cal['header'][i]['OBJECT'] for i in range(len(cal['header']))]
    data = {}

    for key in keys:
        data[key] = {'datetime':[],'optimal_aperture':[],
            'xcen':[],'ycen':[],'fwhm':[],'aspect':[],
            'snrmax':[],'totflux':[],'totflux_aperture':[],
            'chisq':[],'curve_of_growth':[],'secz':[]}

    # keys for optimal aperture dict but not any extra keys
    dkeys = np.array(data.keys())
    dkeys = dkeys[dkeys != 'datetime']

    # go through the files and append name, date as date + UT, and
    # and optimal aperture indices to data

    for i in range(len(files)):
        ra = angcor(cal['header'][i]['RA']).d
        dec = angcor(cal['header'][i]['DEC']).d
        op = tp.optimal_aperture(cal['image'][i],cal['header'][i],[20,30],ra=ra,dec=dec)
        key = cal['header'][i]['OBJECT']
        data[key]['optimal_aperture'] = np.append(data[key]['optimal_aperture'],op['optimal_aperture'])
        data[key]['xcen'] = np.append(data[key]['xcen'],op['xcen'])
        data[key]['ycen'] = np.append(data[key]['ycen'],op['ycen'])
        data[key]['fwhm'] = np.append(data[key]['fwhm'],op['fwhm'])
        data[key]['aspect'] = np.append(data[key]['aspect'],op['aspect'])
        data[key]['snrmax'] = np.append(data[key]['snrmax'],op['snrmax'])
        data[key]['totflux'] = np.append(data[key]['totflux'],op['totflux'])
        data[key]['totflux_aperture'] = np.append(data[key]['totflux_aperture'],op['totflux_aperture'])
        data[key]['chisq'] = np.append(data[key]['chisq'],op['chisq'])
        data[key]['curve_of_growth'] = np.append(data[key]['curve_of_growth'],op['curve_of_growth'])
        data[key]['secz'] = np.append(data[key]['secz'],op['secz'])
        datetime = cal['header'][i]['DATE'] + " " + cal['header'][i]['UT']
        data[key]['datetime'] = np.append(data[key]['datetime'],datetime)

    return data

def pt2(data):

    # format it like reduc fo that it is more robust

    SAmag = np.array([-2.5*np.log10(flux) for flux in data['SA38326']['totflux']])
    GDAmag = np.array([-2.5*np.log10(flux) for flux in data['GD391A']['totflux']])
    GDEmag = np.array([-2.5*np.log10(flux) for flux in data['GD391E']['totflux']])
    KICmag = np.array([-2.5*np.log10(flux) for flux in data['KIC8462852']['totflux']])

    SA = {'mag':SAmag,'secz':data['SA38326']['secz']}
    GDA = {'mag':GDAmag,'secz':data['GD391A']['secz']}
    GDE = {'mag':GDEmag,'secz':data['GD391E']['secz']}
    KIC = {'mag':KICmag,'secz':data['KIC8462852']['secz']}

    SA['fit'] = np.polyfit(SA['secz'],SA['mag'],1)
    GDA['fit'] = np.polyfit(GDA['secz'],GDA['mag'],1)
    GDE['fit'] = np.polyfit(GDE['secz'],GDE['mag'],1)
    KIC['fit'] = np.polyfit(KIC['secz'],KIC['mag'],1)

    SA['b'] = 9.947 - SA['fit'][1]
    GDA['b'] = 12.315 - GDA['fit'][1]
    GDE['b'] = 12.409 - GDE['fit'][1]
    KIC['b'] = 12.01 - KIC['fit'][1]

    data = {'SA':SA,'GDA':GDA,'GDE':GDE,'KIC':KIC}

    return data

def plot_fluxes(data):

    plt.ion()
    plt.figure(1000)
    plt.clf()

    SAflux = data['SA38326']['totflux']
    GDAflux = data['GD391A']['totflux']
    GDEflux = data['GD391E']['totflux']
    KICflux = data['KIC8462852']['totflux']

    SAsecz = data['SA38326']['secz']
    GDAsecz = data['GD391A']['secz']
    GDEsecz = data['GD391E']['secz']
    KICsecz = data['KIC8462852']['secz']

    SAline, = plt.plot(SAsecz,SAflux,'ro',label='SA 38-326')
    GDAline, = plt.plot(GDAsecz,GDAflux,'go',label='GD 391 A')
    GDEline, = plt.plot(GDEsecz,GDEflux,'yo',label='GD 391 E')
    KICline, = plt.plot(KICsecz,KICflux,'bo',label='KIC 8462852')

    plt.ylabel('totflux')
    plt.xlabel('secz')
    ax = plt.gca()
    ax.grid(True)
    plt.title('secz vs total flux on ' + data['GD391A']['datetime'][0][0:8])
    plt.legend(handles=[SAline,GDAline,GDEline,KICline],loc='best')
    plt.show()

    plt.figure(1001)
    plt.clf()

    SAtime = data['SA38326']['datetime']
    GDAtime = data['GD391A']['datetime']
    GDEtime = data['GD391E']['datetime']
    KICtime = data['KIC8462852']['datetime']

    SAtime = [dt.strptime(time,'%d/%m/%y %H:%M:%S') for time in SAtime]
    GDAtime = [dt.strptime(time,'%d/%m/%y %H:%M:%S') for time in GDAtime]
    GDEtime = [dt.strptime(time,'%d/%m/%y %H:%M:%S') for time in GDEtime]
    KICtime = [dt.strptime(time,'%d/%m/%y %H:%M:%S') for time in KICtime]

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
    plt.title('secz vs total flux on ' + data['GD391A']['datetime'][0][0:8])
    plt.legend(handles=[SAline,GDAline,GDEline,KICline],loc='best')
    plt.show()

def plot_mags(data):

    plt.ion()
    plt.figure("secz vs corrected mags")
    plt.clf()

    x = np.linspace(0,2.5,100)

    SAscatter, = plt.plot(data['SA']['secz'],data['SA']['mag']+data['SA']['b'],'ro')
    SAy = [data['SA']['fit'][0] * i + (data['SA']['fit'][1]+data['SA']['b']) for i in x]
    SAmodel, = plt.plot(x,SAy,'r--',label='SA 38-426; mag = '+str(round(data['SA']['fit'][0],3))+'secz'+' + '+str(data['SA']['fit'][1]+data['SA']['b']))

    GDAscatter, = plt.plot(data['GDA']['secz'],data['GDA']['mag']+data['GDA']['b'],'go')
    GDAy = [data['GDA']['fit'][0] * i + (data['GDA']['fit'][1]+data['GDA']['b']) for i in x]
    GDAmodel, = plt.plot(x,GDAy,'g--',label='GD 391 A; mag = '+str(round(data['GDA']['fit'][0],3))+'secz'+' + '+str(data['GDA']['fit'][1]+data['GDA']['b']))

    GDEscatter, = plt.plot(data['GDE']['secz'],data['GDE']['mag']+data['GDE']['b'],'yo')
    GDEy = [data['GDE']['fit'][0] * i + (data['GDE']['fit'][1]+data['GDE']['b']) for i in x]
    GDEmodel, = plt.plot(x,GDEy,'y--',label='GD 391 E; mag = '+str(round(data['GDE']['fit'][0],3))+'secz'+' + '+str(data['GDE']['fit'][1]+data['GDE']['b']))

    KICscatter, = plt.plot(data['KIC']['secz'],data['KIC']['mag']+data['KIC']['b'],'bo')
    KICy = [data['KIC']['fit'][0] * i + (data['KIC']['fit'][1]+data['KIC']['b']) for i in x]
    KICmodel, = plt.plot(x,KICy,'b--',label='KIC 8462852; mag = '+str(round(data['KIC']['fit'][0],3))+'secz'+' + '+str(data['KIC']['fit'][1]+data['KIC']['b']))

    plt.ylabel('corrected mag')
    plt.xlabel('secz')
    plt.title('secz vs corrected mag on 04/12/17')
    ax = plt.gca()
    ax.grid(True)
    plt.legend(handles=[SAmodel,GDAmodel,GDEmodel,KICmodel],loc='best')
    plt.show()
#files, fct = tp.get_files(dir='/Users/nickedwards/python/independent/',prefix='*',tag='solved')
#cal = divide_flat(files)
#data = reduc(files)
#plot_fluxes(data)
#pt2(data)