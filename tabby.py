# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:46:50 2017

@author: nickedwards
"""

# counts per mag = totflux/mag
# mags per airmass = slope from fit
# gain = e- / count
# SA 38-326 mag wrong?
# landolt standard dataframe?


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
import pandas as pd

# values for each star for weighted averages
tot = (1/42.98)+(1/10.22)+(1/22.27)
SAweight = (1/42.98)/tot
GDAweight = (1/10.22)/tot
GDEweight = (1/22.27)/tot

def subtract_bias(files,dir='/Users/nickedwards/python/independent/'):

    bias_image, bias_header = fits.getdata(dir+'master_bias.fits',header=True)

    data = {'image':np.zeros((len(files),2048,2048)), 'header':[]}
    for i in range(len(files)):
        image, header = fits.getdata(files[i],header=True)
        data['image'][i,:,:] = image - bias_image
        data['header'].append(header)

    return data

def divide_flat(files,dir='/Users/nickedwards/python/independent/',fil='V'):

    flat_image, flat_header = fits.getdata(dir+'master_flat_'+fil+'.fits', header= True)

    data = subtract_bias(files)

    for i in range(len(data['image'])): data['image'][i] /= flat_image

    return data

def reduc(files,dir='/Users/nickedwards/python/independent/'):

    cal = divide_flat(files,dir=dir)

    keys = [cal['header'][i]['OBJECT'] for i in range(len(cal['header']))]
    data = {}

    for key in keys:
        data[key] = {'datetime':[],'optimal_aperture':[],
            'xcen':[],'ycen':[],'fwhm':[],'aspect':[],
            'snrmax':[],'totflux':[],'totflux_aperture':[],
            'chisq':[],'curve_of_growth':[],'secz':[],'totfluxerr':[]}

    dkeys = np.array(data.keys())
    dkeys = dkeys[dkeys != 'datetime']

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
        data[key]['totfluxerr'] = np.append(data[key]['totfluxerr'],op['totfluxerr'])
        datetime = cal['header'][i]['DATE'] + " " + cal['header'][i]['UT']
        datetime = dt.strptime(datetime,'%d/%m/%y %H:%M:%S')
        data[key]['datetime'] = np.append(data[key]['datetime'],datetime)

    return data

def pt2(data):

    SAmag = np.array([-2.5*np.log10(flux) for flux in data['SA38326']['totflux']])
    GDAmag = np.array([-2.5*np.log10(flux) for flux in data['GD391A']['totflux']])
    GDEmag = np.array([-2.5*np.log10(flux) for flux in data['GD391E']['totflux']])
    KICmag = np.array([-2.5*np.log10(flux) for flux in data['KIC8462852']['totflux']])

    SA = {'raw mag':SAmag,'secz':data['SA38326']['secz'],
          'totflux':data['SA38326']['totflux'],
          'totfluxerr':data['SA38326']['totfluxerr'],
          'datetime':data['SA38326']['datetime']}
    GDA = {'raw mag':GDAmag,'secz':data['GD391A']['secz'],
           'totflux':data['GD391A']['totflux'],
           'totfluxerr':data['GD391A']['totfluxerr'],
           'datetime':data['GD391A']['datetime']}
    GDE = {'raw mag':GDEmag,'secz':data['GD391E']['secz'],
           'totflux':data['GD391E']['totflux'],
           'totfluxerr':data['GD391E']['totfluxerr'],
           'datetime':data['GD391E']['datetime']}
    KIC = {'raw mag':KICmag,'secz':data['KIC8462852']['secz'],
           'totflux':data['KIC8462852']['totflux'],
           'totfluxerr':data['KIC8462852']['totfluxerr'],
           'datetime':data['KIC8462852']['datetime']}

    SAfit = np.polyfit(SA['secz'],SA['raw mag'],1)
    GDAfit = np.polyfit(GDA['secz'],GDA['raw mag'],1)
    GDEfit = np.polyfit(GDE['secz'],GDE['raw mag'],1)

    SA['m/a'] = np.zeros(len(SA['raw mag']))
    GDA['m/a'] = np.zeros(len(GDA['raw mag']))
    GDE['m/a'] = np.zeros(len(GDE['raw mag']))
    KIC['m/a'] = np.zeros(len(KIC['raw mag']))

    for i in range(len(SA['m/a'])): SA['m/a'][i] = SAfit[0].item()
    for i in range(len(GDA['m/a'])): GDA['m/a'][i] = GDAfit[0].item()
    for i in range(len(GDE['m/a'])): GDE['m/a'][i] = GDEfit[0].item()
    KICma = SAweight*SAfit[0].item()+GDAweight*GDAfit[0].item()+GDEweight*GDEfit[0].item()
    for i in range(len(KIC['m/a'])): KIC['m/a'][i] = KICma

    SA['b'] = np.zeros(len(SA['raw mag']))
    GDA['b'] = np.zeros(len(GDA['raw mag']))
    GDE['b'] = np.zeros(len(GDE['raw mag']))
    KIC['b'] = np.zeros(len(KIC['raw mag']))

    for i in range(len(SA['b'])): SA['b'][i] = 9.947 - SAfit[1]
    for i in range(len(GDA['b'])): GDA['b'][i] = 12.315 - GDAfit[1]
    for i in range(len(GDE['b'])): GDE['b'][i] = 12.409 - GDEfit[1]
    KICb = SAweight*(9.947 - SAfit[1])+GDAweight*(12.315 - GDAfit[1])+GDEweight*(12.409 - GDEfit[1])
    for i in range(len(KIC['b'])): KIC['b'][i] = KICb

    # is this calc right because there should be more counts for SA but fewer
    # mags b/c mags are weird, so the ratio will be off? inversely proportional
    SA['c/m'] = SA['totflux']/9.947
    GDA['c/m'] = GDA['totflux']/12.315
    GDE['c/m'] = GDE['totflux']/12.409
    KIC['c/m'] = KIC['totflux']/12.01
    # different lengths possibly
    # SAweight*(SA['totflux']/9.947)+GDAweight*(GDA['totflux']/12.315)+GDEweight*(GDE['totflux']/12.409)

    # count per mag at 0 airmass

    SA['raw mag'] = SA['b']+SA['raw mag']
    GDA['raw mag'] = GDA['b']+GDA['raw mag']
    GDE['raw mag'] = GDE['b']+GDE['raw mag']
    KIC['raw mag'] = KIC['b']+KIC['raw mag']

    SA['mag0'] = np.zeros(len(SA['raw mag']))
    GDA['mag0'] = np.zeros(len(GDA['raw mag']))
    GDE['mag0'] = np.zeros(len(GDE['raw mag']))
    KIC['mag0'] = np.zeros(len(KIC['raw mag']))

    # b = -ah+k
    # a = m/a
    # h = secz
    # k = raw mag
    for i in range(len(SA['mag0'])): SA['mag0'][i] = SA['raw mag'][i] - SA['m/a'][i]*SA['secz'][i]
    for i in range(len(GDA['mag0'])): GDA['mag0'][i] = GDA['raw mag'][i] - GDA['m/a'][i]*GDA['secz'][i]
    for i in range(len(GDE['mag0'])): GDE['mag0'][i] = GDE['raw mag'][i] - GDE['m/a'][i]*GDE['secz'][i]
    for i in range(len(KIC['mag0'])): KIC['mag0'][i] = KIC['raw mag'][i] - KIC['m/a'][i]*KIC['secz'][i]

    data = {'SA38326':SA,'GD391A':GDA,'GD391E':GDE,'KIC8462852':KIC}

    return data

def standards(data):

    SAmag = np.array([-2.5*np.log10(flux) for flux in data['SA38326']['totflux']])
    GDAmag = np.array([-2.5*np.log10(flux) for flux in data['GD391A']['totflux']])
    GDEmag = np.array([-2.5*np.log10(flux) for flux in data['GD391E']['totflux']])

    SA = {'raw mag':SAmag,'secz':data['SA38326']['secz'],
          'totflux':data['SA38326']['totflux'],
          'totfluxerr':data['SA38326']['totfluxerr'],
          'datetime':data['SA38326']['datetime']}
    GDA = {'raw mag':GDAmag,'secz':data['GD391A']['secz'],
           'totflux':data['GD391A']['totflux'],
           'totfluxerr':data['GD391A']['totfluxerr'],
           'datetime':data['GD391A']['datetime']}
    GDE = {'raw mag':GDEmag,'secz':data['GD391E']['secz'],
           'totflux':data['GD391E']['totflux'],
           'totfluxerr':data['GD391E']['totfluxerr'],
           'datetime':data['GD391E']['datetime']}

    SAfit = np.polyfit(SA['secz'],SA['raw mag'],1)
    GDAfit = np.polyfit(GDA['secz'],GDA['raw mag'],1)
    GDEfit = np.polyfit(GDE['secz'],GDE['raw mag'],1)

    SA['m/a'] = SAfit[0].item()
    GDA['m/a'] = GDAfit[0].item()
    GDE['m/a'] = GDEfit[0].item()

    SA['b'] = 9.947 - SAfit[1]
    GDA['b'] = 12.315 - GDAfit[1]
    GDE['b'] = 12.409 - GDEfit[1]

    # b = -ah+k
    # a = m/a
    # h = secz
    # k = raw mag + b

    """
    use these lines in your script.  I use a package called astroquery.simbad.
    from astroquery.simbad import Simbad
    Simbad.add_votable_fields('sptype')
    Simbad.add_votable_fields('measurements')
    then to query by starname use Simbad.query_object(starname).to_pandas()
    """
    SA['mag'] = SA['raw mag'] + 9.947 - SAfit[1]
    GDA['mag'] = GDA['raw mag'] + 12.315 - GDAfit[1]
    GDE['mag'] = GDE['raw mag'] + 12.409 - GDEfit[1]

    SAcmfit = np.polyfit(SA['secz'],SA['totflux'],1)
    GDAcmfit = np.polyfit(GDA['secz'],GDA['totflux'],1)
    GDEcmfit = np.polyfit(GDE['secz'],GDE['totflux'],1)

    SA['c0'] = SAcmfit[1].item()*(10**(.4*9.947))
    GDA['c0'] = GDAcmfit[1].item()*(10**(.4*12.315))
    GDE['c0'] = GDEcmfit[1].item()*(10**(.4*12.409))

    SA['mag0'] = SA['mag'] - SA['m/a']*SA['secz']
    GDA['mag0'] = GDA['mag'] - GDA['m/a']*GDA['secz']
    GDE['mag0'] = GDE['mag'] - GDE['m/a']*GDE['secz']

    data = {'SA38326':SA,'GD391A':GDA,'GD391E':GDE}
    return data

def source(data, standard):

    KIC = {'secz':data['KIC8462852']['secz'],
           'totflux':data['KIC8462852']['totflux'],
           'totfluxerr':data['KIC8462852']['totfluxerr'],
           'datetime':data['KIC8462852']['datetime']}

    KIC['raw mag'] = np.array([-2.5*np.log10(flux) for flux in KIC['totflux']])

    b = SAweight*standard['SA38326']['b']+GDAweight*standard['GD391A']['b']+GDEweight*standard['GD391E']['b']
    ma = SAweight*standard['SA38326']['m/a']+GDAweight*standard['GD391A']['m/a']+GDEweight*standard['GD391E']['m/a']

    KIC['mag0'] = KIC['raw mag']+b - ma*KIC['secz']
    #brighter than expected probably because SA is weird
    return KIC

def save(data,old=True):
    frame = pd.DataFrame.from_dict(data)
    if old:
        old = pd.read_csv('/Users/nickedwards/python/independent/tabby.csv')
        # check to see if you reassign dataframes if index follow over
        # old = old[old['datetime']]

    return frame

def plot_fluxes(data):

    plt.ion()
    for key in data.keys():
        plt.figure('secz vs total flux for ' + key)
        plt.clf()

        flux = data[key]['totflux']

        secz = data[key]['secz']

        line, = plt.plot(secz,flux,'bo',label=key)

        plt.ylabel('totflux')
        plt.xlabel('secz')
        ax = plt.gca()
        ax.grid(True)
        plt.title('secz vs total flux for ' + key)
        plt.legend(handles=[line],loc='best')
        plt.show()

        plt.figure('secz vs time for ' + key)
        plt.clf()

        time = data[key]['datetime']

        line, = plt.plot_date(time,flux,'bo',label=key)

        ax = plt.gca()
        xfmt = DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(xfmt)
        ax.grid(True)
        plt.ylabel('totflux')
        plt.xlabel('time (UT)')
        plt.title('secz vs time for ' + key)
        plt.legend(handles=[line],loc='best')
        plt.show()

def plot_mags(data):

    plt.ion()
    for key in data.keys():
        plt.figure("secz vs corrected mags for " +str(key))
        plt.clf()

        x = np.linspace(0,2.5,100)

        scatter, = plt.plot(data[key]['secz'],data[key]['raw mag'],'bo')
        y = [data[key]['m/a'][0] * i + data[key]['mag0'][0] for i in x]
        model, = plt.plot(x,y,'b--',label=key+'; mag = '+str(round(data[key]['m/a'][0],3))+'secz'+' + '+str(data[key]['mag0'][0]))

        plt.ylabel('corrected mag')
        plt.xlabel('secz')
        plt.title('secz vs corrected mag for ' + key)
        ax = plt.gca()
        ax.grid(True)
        plt.legend(handles=[model],loc='best')
        plt.show()
"""
files, fct = tp.get_files(dir='/Users/nickedwards/python/independent/',prefix='*',tag='solved')
cal = divide_flat(files)
data = reduc(files)
stuff = pt2(data)
standard = standards(data)
csv = save(stuff['KIC8462852'])
"""
"""
import matplotlib; matplotlib.use('Agg')
import thacherphot as tp
files, fct = tp.get_files(dir='/home/student/nedwards/2017Apr25/',suffix='V.fts')
tp.do_astrometry(files)
"""