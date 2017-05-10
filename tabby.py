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

    data = subtract_bias(files,dir=dir)

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

    SAcmfit = np.polyfit(SA['secz'],SA['totflux']/10,1)
    GDAcmfit = np.polyfit(GDA['secz'],GDA['totflux']/90,1)
    GDEcmfit = np.polyfit(GDE['secz'],GDE['totflux']/100,1)

    SA['c0'] = SAcmfit[1].item()*(10**(.4*9.947))
    GDA['c0'] = GDAcmfit[1].item()*(10**(.4*12.315))
    GDE['c0'] = GDEcmfit[1].item()*(10**(.4*12.409))

    SA['Mzp'] = -2.5*np.log10(1/SA['c0'])
    GDA['Mzp'] = -2.5*np.log10(1/GDA['c0'])
    GDE['Mzp'] = -2.5*np.log10(1/GDE['c0'])

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

    mzpave = (standard['SA38326']['Mzp']+standard['GD391A']['Mzp']+standard['GD391E']['Mzp'])/3
    KIC['raw mag ave'] = mzpave-2.5*np.log10(KIC['totflux']/60)
    maave = (standard['SA38326']['m/a']+standard['GD391A']['m/a']+standard['GD391E']['m/a'])/3

    KIC['mag0 ave'] = KIC['raw mag ave'] - maave*KIC['secz']
    return KIC

def save(data):
    tabby = pd.DataFrame.from_dict(data)
    tabby.to_csv('/home/student/nedwards/independent/')
    return frame

def finish():
    files, fct = tp.get_files(dir='/home/student/nedwards/2017A*/',suffix='solved.fits')
    data = reduc(files,dir='/home/student/nedwards/independent/')
    standard = standards(data)
    tabby = source(stuff,standard)
    save(tabby)
    return

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