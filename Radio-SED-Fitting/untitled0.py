#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 19:42:28 2020

@author: pallavipatil
"""

# SED Analysis Plots:

import sys
sys.path.insert(1, '/Users/pallavipatil/Desktop/VLA/Radio-SED-Fitting/')
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table,join
import os
import seaborn as sns
import matplotlib.image as mpimg
import Radio_Models_func
import file_prep_radio_fitting 
from matplotlib.patches import Rectangle
import fitting_mcsim_func as mcsim_fit
import json
import pandas as pd


loadt = file_prep_radio_fitting.LoadTables()
fprep = file_prep_radio_fitting.FilePrepUtils()
rmfit = Radio_Models_func.RadioModelFit()
atscat, atgcat, vla_ax, vla_bx = loadt.get_tabs()

wise_phot = Table.read('WISE_phot_comb_manual_fluxes.csv')
wise_phot.add_index('Source_name')

vla_ax.add_index('Source_name')
ind = vla_ax.loc_indices['J0643-27']
if ind>0:
    vla_ax.remove_row(ind)
    

if 'WISEname' not in vla_ax.colnames:
    vla_ax = join(vla_ax, wise_phot['Source_name', 'WISEname'], keys='Source_name', 
              join_type='left')

if 'WISEname' not in vla_bx.colnames:
    vla_bx = join(vla_bx, wise_phot['Source_name', 'WISEname'], keys='Source_name', 
              join_type='left')

vla_ax_grp = vla_ax.group_by('WISEname')
vla_bx_grp = vla_bx.group_by('WISEname')
#####################3
# File with spectral class info:


spdf = pd.read_excel('New_Analysis/Spectral_shape_classf_master.xlsx')
sp_class = Table.from_pandas(spdf)
sp_class.add_index('Source_name')

final_unq_jvla = Table.read('New_Analysis/final_unq_JVLASP_vf_morph_regmod.csv')
final_unq_jvla.add_index('Source_name')

fit_results_tab = Table.read('New_Analysis/fit_sed_results_tab.csv')
fit_results_tab.add_index('source')

#final_tab_jvla = Table.read('New_Analysis/final_class_JVLASP_vf_morph_regmod.csv')
final_tab_jvla = Table.read('New_Analysis/final_tab_new_jvla.csv')
final_tab_jvla.add_index('Source_name')


fjvla_AX = Table.read('../VLA-NewWork/JMFIT_CASA_AX.csv')
fjvla_BX = Table.read('../VLA-NewWork/JMFIT_CASA_BX.csv')

fjvla_AX.add_index('Source_name')
fjvla_BX.add_index('Source_name')

fp = open('cpl_fits.json')
par_fits = json.load(fp)

atgcat.add_index('WISEname')
atscat.add_index('WISEname')
sp_class.add_index('Source_name')

survey_info = dict.fromkeys(['NVSS', 'FIRST', 'SUMSS', 'TGSS', 'VLSSr', 'WENSS',
                          'WISH', 'GB6', 'VLITE', 'VLASS', 'TEXAS',  'GLEAM', 'LOTSS', 
                          'AX', 'BX'])
freq_survey = [1.4, 1.4, 0.843, 0.150, 0.074, 0.325, 0.352, 4.85, 0.338, 3.0, 0.320, 
                           0.200, 0.150, 10, 10]

for key, freq in zip(survey_info.keys(), freq_survey):
    fdict = {'Freq':freq}
    survey_info[key] = fdict
    
# Radio Color-Color Plot: it plots alpha-high and alpha-low:
# If TGSS is available: alpha-low will be calculated between NVSS and TGSS 
# and alpha-high: VLA-X and NVSS/VLASS:

alhigh = []
allow = []
e_alhigh = []
e_allow = []
alow_lim = []

flux_dict = dict.fromkeys(wise_phot['Source_name'])

for wise_row in wise_phot:
    source = wise_row['Source_name']
    wisen = wise_row['WISEname']
    sp_row = sp_class.loc[source]
    glmcat = atgcat.loc[wisen]
    uplim = False
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
         
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
            OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, wise_row, jvla_AX,jvla_BX)
    freq_arr, flux_arr, eflux_arr, labels = mcsim_fit.modify_data(source,freq_arr,
                            flux_arr, eflux_arr, labels,sp_row, wise_row )
    pardict = dict.fromkeys(labels)
    for nu, f, ef, label in zip(freq_arr, flux_arr, eflux_arr, labels):
        pardict[label] = [nu, f, ef]
    flux_dict[source] = pardict

    hcat = 'BX'
    lcat = 'TGSS'
    if 'BX' not in labels:
        hcat = 'AX'
    if 'TGSS' not in labels:
        cind = np.argmin(freq_arr[freq_arr>0])
        lcat = labels[cind]

    mycat = [lcat, 'NVSS', hcat]
    
    # check for presence of limit on lower frequency data
    lcol = lcat+'_Limit'
    if lcol in sp_row.colnames:
        if sp_row[lcol] == 'U':
            uplim = True
    # calculation of alpha_high
    s1 =         
            
    alow_lim.append(uplim)


