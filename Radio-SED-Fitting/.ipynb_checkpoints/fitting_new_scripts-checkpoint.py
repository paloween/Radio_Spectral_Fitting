#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:37:12 2019

@author: pallavipatil
"""

import numpy as np
import pandas as pd
from astropy.table import Table, join, Column
from astropy.io import fits
import os
import errors_mcsim as emcsim
import sys
sys.path.insert(1,'../Radio-SED-Fitting')
import Radio_Models_func
import file_prep_radio_fitting 
import fitting_mcsim_func as mcsim_fit

loadt = file_prep_radio_fitting.LoadTables()
fprep = file_prep_radio_fitting.FilePrepUtils()
rmfit = Radio_Models_func.RadioModelFit()

##Accessing all different data files and creating astropy tables.

atscat, atgcat, vla_ax, vla_bx = loadt.get_tabs()

wise_phot  = Table.read('WISE_phot_comb_manual_fluxes.csv')
wise_phot.add_index('WISEname')
wise_phot.add_index('Source_name')

if 'WISEname' not in vla_ax.colnames:
    vla_ax = join(vla_ax, wise_phot['Source_name', 'WISEname'], keys='Source_name', 
              join_type='left')

if 'WISEname' not in vla_bx.colnames:
    vla_bx = join(vla_bx, wise_phot['Source_name', 'WISEname'], keys='Source_name', 
              join_type='left')
    
vla_ax.add_index('Source_name')
ind = vla_ax.loc_indices['J0643-27']
vla_ax.remove_row(ind)   
vla_ax.remove_indices('Source_name')

vla_ax_grp = vla_ax.group_by('WISEname')
vla_bx_grp = vla_bx.group_by('WISEname')
    
sp_class = loadt.get_sp_info()
sp_class_sources = list(sp_class[sp_class['Obs_Flag']!='NO']['Source_name'])
sp_class.add_index('Source_name')

spdf = pd.read_excel('New_Analysis/Spectral_shape_classf_master.xlsx')
sp_class = Table.from_pandas(spdf)
sp_class.add_index('Source_name')



# fitTab = ascii.read('GuessPar_Radiofits.csv', format='basic', delimiter=',')
# fitTab.add_index('Name')

fit_results_ext_tab = Table(names = ('source', 'region',  'PL_s0', 'PL_al', 'CPL_s0', 
                                'CPL_al','CPL_nup',  'EFFA_s0', 'EFFA_al', 'EFFA_nup',
                                 'IFFA_s0', 'IFFA_al', 'IFFA_nup','SSA_s0', 
                                 'SSA_al', 'SSA_nup', 'perr_pl_s0', 'perr_pl_al',
                                 'perr_cpl_s0','perr_cpl_al', 'perr_cpl_nup', 
                                 'perr_effa_s0','perr_effa_al', 'perr_effa_nup',
                                 'perr_iffa_s0','perr_iffa_al', 'perr_iffa_nup',
                                 'perr_ssa_s0','perr_ssa_al', 'perr_ssa_nup',
                                 'rchisq_pl', 'rchisq_cpl','rchisq_effa',
                                 'rchisq_iffa','rchisq_ssa'), 
                        dtype = ('S10', 'S10', 'f8', 'f8', 'f8', 'f8', 'f8',
                                 'f8', 'f8', 'f8','f8', 'f8', 'f8','f8', 'f8', 
                                 'f8','f8','f8','f8','f8','f8','f8','f8','f8',
                                 'f8','f8','f8','f8','f8','f8','f8','f8','f8',
                                 'f8','f8'))
    
cpl_fits = dict.fromkeys(wise_phot['Source_name'])
for key, value in cpl_fits.items():
    cpl_fits[key] = dict.fromkeys(['s0_cpl', 'al_cpl', 'q_cpl', 'speak', 'nupeak' ])
    

mcsim_flag = True
wise_phot.sort('WISEname')

for wise_row in wise_phot:
    best_fits = []
    chi_sqs = []
    par_err = []
    mcsim = []
    wisen = wise_row['WISEname']
    scat = wise_row
    glmcat = atgcat.loc[wisen]
    sname = wise_row['Source_name']
    sp_row = sp_class.loc[sname]
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
    OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,
                                           jvla_BX)
    freq_arr, flux_arr, eflux_arr, labels = mcsim_fit.modify_data(sname,freq_arr, flux_arr,
                                eflux_arr, labels,sp_row, scat  )
    
    if 'BX' in labels and np.any(jvla_BX['snr']>50):
        freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, 
                                        flux_arr, eflux_arr, alpha_BX)
    if 'AX' in labels and np.any(jvla_AX['snr']>50):
        freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, 
                                        flux_arr, eflux_arr, alpha_AX)

   
    # Guess parameters:
    #guess_pars = fitTab.loc[wisen]
    models = ['PL', 'CPL', 'EFFA', 'IFFA', 'SSA']
    # Perform radio fitting:
    myrow = [sname, 'reg1' ]
    for model in models:
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
        
        fit_res, perr_res, chi_res = rmfit.chi_sq_fit(freq_arr, flux_arr,
                     eflux_arr, guess_pars, model)
        best_fits.append(fit_res)
        par_err.extend(perr_res)
        chi_sqs.append(chi_res)
        myrow.extend(list(fit_res))
    if mcsim_flag:
        pars, par_fits, tcarlo  = mcsim_fit.mcsim_fitting(wise_row, 'CPL', atgcat, vla_ax_grp, 
                                                            vla_bx_grp, 1000)
        for row, key in zip(par_fits, ['s0_cpl', 'al_cpl', 'q_cpl', 'speak', 'nupeak' ]):
            mcsim.append(row)
            pardict = {} 
            for k, v in zip(pars, row):
                pardict[k] = v
            cpl_fits[sname][key] = pardict            
    #print(sname, cpl_fits[sname]['s0_cpl']['Mean'], cpl_fits[sname]['al_cpl']['Mean'],cpl_fits[sname]['q_cpl']['Mean'] )  
    #print(sname, fit_res )        
    myrow.extend(par_err)
    myrow.extend(chi_sqs)
    fit_results_ext_tab.add_row(myrow)
    print(sname)
 
    
# Cleaning data products:
        
      
wise_phot.add_column(Column(data = [np.nan]*len(wise_phot), name = 'TGSS_Limit',
                            dtype='f8'))
tgss_limits = sp_class[sp_class['TGSS_Limit']=='U']['Source_name']
tgss_dir = '/Users/pallavipatil/Desktop/VLA/Radio_Fits_v2/Astroquery/TGSS/'


for name in tgss_limits:
    tgss_file = name+'_TGSS.fits'
    if os.path.exists(tgss_dir+tgss_file):
        hdu = fits.open(tgss_dir+tgss_file)
        data = hdu[0].data
        header = hdu[0].header
        tgss_clip  =np.nanpercentile(data, 99)
        tgss_rms = np.nanstd(data[data<tgss_clip])
        print(name, tgss_rms)
        wise_ind = wise_phot.loc_indices[name]
        wise_phot['TGSS_Limit'][wise_ind] = tgss_rms*1000
               

wise_phot.add_column(Column(data = [np.nan]*len(wise_phot), name = 'VLSSr_Limit',
                            dtype='f8'))

vlssr_limits = sp_class[sp_class['VLSSr_Limit']=='U']['Source_name']
vlssr_dir = '/Users/pallavipatil/Desktop/VLA/Radio_Fits_v2/Astroquery/VLSSr/'

for name in vlssr_limits:
    vlssr_file = name+'_VLSSr.fits'
    if os.path.exists(vlssr_dir+vlssr_file):
        hdu = fits.open(vlssr_dir+vlssr_file)
        data = hdu[0].data
        header = hdu[0].header
        vlssr_clip  =np.nanpercentile(data, 99)
        vlssr_rms = np.nanstd(data[data<vlssr_clip])
        print(name, vlssr_rms)
        wise_ind = wise_phot.loc_indices['Source_name',name]
        wise_phot['VLSSr_Limit'][wise_ind] = vlssr_rms*1000

wise_phot.add_column(Column(data = [np.nan]*len(wise_phot), name = 'GLEAM_Limit',
                            dtype='f8'))

gleam_limits = sp_class[sp_class['GLEAM_Limit']=='U']['Source_name']
gleam_dir = '/Users/pallavipatil/Desktop/VLA/Radio_Fits_v2/Astroquery/GLEAM/'

for name in gleam_limits:
    gleam_file = name+'_GLEAM_170-231_10deg.fits'
    if os.path.exists(gleam_dir+gleam_file):
        hdu = fits.open(gleam_dir+gleam_file)
        data = hdu[0].data
        header = hdu[0].header
        gleam_clip  =np.nanpercentile(data, 99)
        gleam_rms = np.nanstd(data[data<gleam_clip])
        print(name, gleam_rms)
        wise_ind = wise_phot.loc_indices['Source_name',name]
        wise_phot['GLEAM_Limit'][wise_ind] = gleam_rms*1000

wise_phot.add_column(Column(data = [np.nan]*len(wise_phot), name = 'WENSS_Limit',
                            dtype='f8'))


wenss_limits = sp_class[sp_class['WENSS_Limit']=='U']['Source_name']
wenss_dir = '/Users/pallavipatil/Desktop/VLA/Radio_Fits_v2/Astroquery/WENSS/'

for name in wenss_limits:
    wenss_file = name+'_WENSS.fits'
    if os.path.exists(wenss_dir+wenss_file):
        hdu = fits.open(wenss_dir+wenss_file)
        data = hdu[0].data
        header = hdu[0].header
        wenss_clip  =np.nanpercentile(data, 99)
        wenss_rms = np.nanstd(data[data<wenss_clip])
        print(name, wenss_rms)
        wise_ind = wise_phot.loc_indices['Source_name',name]
        wise_phot['WENSS_Limit'][wise_ind] = wenss_rms*1000
