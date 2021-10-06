#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:10:22 2020

@author: pallavipatil
"""

import numpy as np
import errors_mcsim as emcsim
import sys
sys.path.insert(1,'../Radio-SED-Fitting')
import Radio_Models_func
import file_prep_radio_fitting 
import pandas as pd
from astropy.table import Table
import os 
from astropy.io import fits

loadt = file_prep_radio_fitting.LoadTables()
fprep = file_prep_radio_fitting.FilePrepUtils()
rmfit = Radio_Models_func.RadioModelFit()
atscat, atgcat, vla_ax, vla_bx = loadt.get_tabs()

spdf = pd.read_excel('../Radio_Fits_v2/New_Analysis/Spectral_shape_classf_master.xlsx')
sp_class = Table.from_pandas(spdf)
sp_class.add_index('Source_name')
    

def modify_data(source, freq_arr, flux_arr, eflux_arr, labels, sp_row, wise_row):
    vlass_flag = sp_row['VLASS_Limit']
    tgss_flag = sp_row['TGSS_Limit']
    vlssr_flag = sp_row['VLSSr_Limit']
    wenss_flag = sp_row['WENSS_Limit']
    gleam_flag = sp_row['GLEAM_Limit']
    sumss_flag = sp_row['SUMSS_Limit']
    first_flag = sp_row['FIRST_Limit']
    freq_arr = list(freq_arr)
    flux_arr = list(flux_arr)
    eflux_arr = list(eflux_arr)
    ul_limit = 7 # in sigma
    if freq_arr.count(10) ==2:
        ind = labels.index('AX')
        del labels[ind]
        del freq_arr[ind]
        del flux_arr[ind]
        del eflux_arr[ind]
    if vlass_flag == 'F' and 'VLASS' in labels:
        ind = labels.index('VLASS')
        del labels[ind]
        del freq_arr[ind]
        del flux_arr[ind]
        del eflux_arr[ind]
        
    if vlass_flag == 'M':
        try:
            ind = labels.index('VLASS')
            flux_arr[ind] = wise_row['FVLASS_man']*1000
            eflux_arr[ind] = wise_row['EVLASS_man']*1000
        except:
            labels.append('VLASS')
            flux_arr.append(wise_row['FVLASS_man']*1000)
            eflux_arr.append(wise_row['EVLASS_man']*1000)
            freq_arr.append(3.0)
    if vlass_flag == 'I':
        ind = labels.index('VLASS')
        flux_arr[ind] = wise_row['Total_flux']*1000
        eflux_arr[ind] = wise_row['E_Total_flux']*1000
        
    if 'TGSS' in labels:
        ind = labels.index('TGSS')
        flux = flux_arr[ind]
        eflux = eflux_arr[ind]
        if eflux<0.2*flux:
            eflux_arr[ind] = flux*0.2
    if tgss_flag == 'U':
        labels.append('TGSS')
        flux_level = ul_limit/2
        flux_arr.append(flux_level*wise_row['TGSS_Limit'])
        eflux_arr.append(flux_level*wise_row['TGSS_Limit'])
        freq_arr.append(0.150)
    if vlssr_flag == 'U':
        noise = wise_row['VLSSr_Limit']
        if noise>0:
            labels.append('VLSSr')
            flux_level = ul_limit/2
            flux_arr.append(flux_level*wise_row['VLSSr_Limit'])
            eflux_arr.append(flux_level*wise_row['VLSSr_Limit'])
            freq_arr.append(0.074)
    if wenss_flag == 'U':
        noise = wise_row['WENSS_Limit']
        if noise>0:
            labels.append('WENSS')
            flux_level = ul_limit/2
            flux_arr.append(flux_level*wise_row['WENSS_Limit'])
            eflux_arr.append(flux_level*wise_row['WENSS_Limit'])
            freq_arr.append(0.325)
    if gleam_flag == 'U':
        noise = wise_row['GLEAM_Limit']
        if noise>0:
            labels.append('GLEAM')
            flux_level = 5/2
            flux_arr.append(flux_level*wise_row['GLEAM_Limit'])
            eflux_arr.append(flux_level*wise_row['GLEAM_Limit'])
            freq_arr.append(0.200)
    if sumss_flag == 'R' and 'SUMSS' in labels:
        ind = labels.index('SUMSS')
        del labels[ind]
        del freq_arr[ind]
        del flux_arr[ind]
        del eflux_arr[ind]
    if first_flag == 'R' and 'FIRST' in labels:
        ind = labels.index('FIRST')
        del labels[ind]
        del freq_arr[ind]
        del flux_arr[ind]
        del eflux_arr[ind]
    if tgss_flag == 'T' and 'TGSS' not in labels:
        flux, eflux = get_faint_detection(source, 'TGSS',  10)
        if flux>0:
            labels.append('TGSS')
            flux_arr.append(flux)
            eflux_arr.append(eflux)
            freq_arr.append(0.15)
    if vlssr_flag == 'F' and 'VLSSr' not in labels:
        flux, eflux = get_faint_detection(source, 'VLSSr', 10)
        if flux>0:
            labels.append('VLSSr')
            flux_arr.append(flux)
            eflux_arr.append(eflux)
            freq_arr.append(0.074)
            
    if gleam_flag == 'F' and 'GLEAM' not in labels:
        flux, eflux = get_faint_detection(source, 'GLEAM',  10)
        if flux>0:
            labels.append('GLEAM')
            flux_arr.append(flux)
            eflux_arr.append(eflux)
            freq_arr.append(0.200)
    if sumss_flag == 'F' and 'SUMSS' not in labels:
        flux, eflux = get_faint_detection(source, 'SUMSS', 10)
        if flux>0:
            labels.append('SUMSS')
            flux_arr.append(flux)
            eflux_arr.append(eflux)
            freq_arr.append(0.843)
    if wenss_flag == 'F' and 'WENSS' not in labels:
        flux, eflux = get_faint_detection(source, 'WENSS', 10)
        if flux>0:
            labels.append('WENSS')
            flux_arr.append(flux)
            eflux_arr.append(eflux)
            freq_arr.append(0.325)
    if gleam_flag == 'C' and 'GLEAM' in labels:
        ind = labels.index('GLEAM')
        del labels[ind]
        del freq_arr[ind]
        del flux_arr[ind]
        del eflux_arr[ind]
    if wenss_flag == 'C' and 'WENSS' in labels:
        ind = labels.index('WENSS')
        del labels[ind]
        del freq_arr[ind]
        del flux_arr[ind]
        del eflux_arr[ind]

        
            
    freq_arr = np.array(freq_arr)
    flux_arr = np.array(flux_arr)
    eflux_arr = np.array(eflux_arr)
    return freq_arr, flux_arr, eflux_arr, labels


def get_faint_detection(source, survey, reg_pix = 10):
    fits_dir = '/Users/pallavipatil/Desktop/VLA/Radio_Fits_v2/Astroquery/'
    fitsf = fits_dir+survey+'/'+source+'_'+survey+'.fits'
    if survey == 'GLEAM':
        fitsf = fits_dir+survey+'/'+source+'_'+survey+'_170-231_10deg.fits'
    if os.path.exists(fitsf):
        hdu = fits.open(fitsf)
        data = hdu[0].data
        head = hdu[0].header
        ylen = head['NAXIS2']
        xlen = head['NAXIS1']
        xcen = int(xlen/2)
        ycen = int(ylen/2)
        xlow = int(xcen-reg_pix/2)
        ylow = int(ycen-reg_pix/2)
        xupp = int(xcen+reg_pix/2)
        yupp = int(ycen+reg_pix/2)
        dslice = data[xlow:xupp, ylow:yupp]
        maxflux = np.max(dslice)
        # calculate rms:
        data_clip = abs(np.nanpercentile(data, 0.001))
        data_rms = np.nanstd(data[data<data_clip])
        return maxflux*1000, data_rms*1000
    else:
        return -999, -999

        
def mcsim_fitting(row, model, atgcat, vla_ax_grp, vla_bx_grp, ntrial, useInBand=True):
    wisen = row['WISEname']
    scat = row
    glmcat = atgcat.loc[wisen]
    sname = row['Source_name']
    sp_row = sp_class.loc[sname]
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
    OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,
                                           jvla_BX)
    
    
    freq_arr, flux_arr, eflux_arr, labels = modify_data(sname,freq_arr, flux_arr,
                                eflux_arr, labels,sp_row, scat )
    
    if useInBand:
        if 'BX' in labels and np.any(jvla_BX['snr']>50):
            freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, 
                                            flux_arr, eflux_arr, alpha_BX)
        if 'AX' in labels and np.any(jvla_AX['snr']>50):
            freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, 
                                            flux_arr, eflux_arr, alpha_AX)

    
   
    # Guess parameters:
    #guess_pars = rmfit.estimate_guess_pars(freq_arr, flux_arr, model)
    #if np.isnan(guess_pars).any():
    s0 = np.max(flux_arr)
    nu_t = np.min(freq_arr[freq_arr > 0])
    alpha = -0.7
    guess_cpl = [s0, alpha, nu_t]
    guess_pl = [s0, alpha]
    if model == 'PL':
        guess_pars = guess_pl
    elif model == 'CPL':
        guess_pars = [s0, alpha, 0.0]
    elif model == 'GCV':
        guess_pars = [s0, -alpha, alpha, nu_t]
    else:
        guess_pars = guess_cpl

    # Perform radio fitting:
    pars, par_fits, tcarlo =  emcsim.bootstrap_mcsim(freq_arr, flux_arr, 
                                        eflux_arr,  guess_pars, model, ntrial)
    
    return pars, par_fits, tcarlo
