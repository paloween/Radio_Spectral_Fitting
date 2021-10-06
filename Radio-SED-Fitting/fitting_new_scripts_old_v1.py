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
import json
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import dquants

def edit_fits_header(fitsfile):
    hdu = fits.open(fitsfile)
    hdr = hdu[0].header
    mydata= hdu[0].data
    survey = os.path.basename(fitsfile).split('_')[1].split('.')[0]
    keywords = ['CTYPE3','CRVAL3','CDELT3', 'CRPIX3', 'CROTA3','CUNIT3',
                'CTYPE4','CRVAL4','CDELT4', 'CRPIX4', 'CROTA4','CUNIT4',
                'PC3_1', 'PC4_1', 'PC3_2','PC4_2', 'PC1_3','PC2_3','PC3_3',
                'PC4_3','PC1_4','PC2_4','PC3_4','PC4_4']
    for key in keywords:
        hdr.remove(key, ignore_missing=True)
    if survey =='NVSS' or 'VLASS':
        hdr.remove('NAXIS3',ignore_missing=True )
        hdr.remove('NAXIS4',ignore_missing=True )
        hdr['NAXIS'] = 2
    if len(mydata.shape)>2:
        mydata = mydata[0,0,:,:]
    new_fits = fits.PrimaryHDU(data=mydata, header=hdr)
    hdu.close()
    return new_fits
    
    
def get_faint_detection(source, survey, fitsf, pos, reg_pix = 10):
    if os.path.exists(fitsf):
        hdu = edit_fits_header(fitsf)
        data = hdu.data
        head = hdu.header
        wcs = WCS(head)
        xcen, ycen = wcs.world_to_pixel(pos)
        xlow = int(xcen-reg_pix/2)
        ylow = int(ycen-reg_pix/2)
        xupp = int(xcen+reg_pix/2)
        yupp = int(ycen+reg_pix/2)
        dslice = data[xlow:xupp, ylow:yupp]
        maxflux = np.max(dslice)
        # calculate rms:
        #data_clip = abs(np.nanpercentile(data, 0.001))
        #data_rms = np.nanstd(data[data<data_clip])
        rms = head['ACTNOISE']
        return maxflux*1000, rms*1000
    else:
        return -999, -999
    
    
def get_VCSS_data(source, wise_row, sp_row, labels, freq_arr, flux_arr, eflux_arr):
    vlite_sources = {'J1238+52':17,  'J0823-06':16.2, 'J0642-27':20 }

    if 'VLITE' in labels:
        my_ind = labels.index('VLITE')
        if eflux_arr[my_ind]/flux_arr[my_ind] <0.2:
            eflux_arr[my_ind] = np.sqrt((0.2*flux_arr[my_ind])**2+eflux_arr[my_ind]**2)
            
    else:
        if sp_row['VLITE_Limit'] == 'F':
            mypos = SkyCoord(ra=wise_row['wise_ra_1'], dec = wise_row['wise_dec_1'], 
                             unit = (u.deg, u.deg), frame='fk5')
            fitsf = './Astroquery/VCSS/cutouts/'+source+'_VCSS1_cutout_v1.fits'
            flux, rms = get_faint_detection(source, 'VLITE', fitsf, mypos, 
                                            reg_pix = 10)
            eflux = np.sqrt(0.2*flux**2+rms**2)
            freq_arr = np.append(freq_arr, 0.338)
            flux_arr = np.append(flux_arr, flux)
            eflux_arr = np.append(eflux_arr, eflux)
            labels.append('VLITE')
        elif wise_row['VLITE_Limit']>0:
            ul_limit = 10/2
            flux = ul_limit*wise_row['VLITE_Limit']
            eflux = ul_limit*wise_row['VLITE_Limit']           
            freq_arr = np.append(freq_arr, 0.338)
            flux_arr = np.append(flux_arr, flux)
            eflux_arr = np.append(eflux_arr, eflux)
            labels.append('VLITE')
 
    #add VLITE data:
    if source in list(vlite_sources.keys()):
        freq_arr = np.append(freq_arr, 0.338)
        flux_arr = np.append(flux_arr, vlite_sources[source])
        eflux_arr = np.append(eflux_arr, vlite_sources[source]*0.2)
        labels.append('VLITE')
    return freq_arr, flux_arr, eflux_arr, labels

def VCSS_data_cleanup(source, wise_row, sp_row, labels, freq_arr, flux_arr, eflux_arr):
    
    if 'VLITE' in labels:
        my_ind = labels.index('VLITE')
        fitsf = './Astroquery/VCSS/cutouts/'+source+'_VCSS1_cutout_v1.fits'
        header = fits.open(fitsf)[0].header
        snr = flux_arr[my_ind]/header['ACTNOISE']/1000
        if snr<10 or source == 'J0526-32':
            freq_arr = np.delete(freq_arr, my_ind)
            flux_arr = np.delete(flux_arr, my_ind)
            eflux_arr = np.delete(eflux_arr, my_ind)
            labels.remove('VLITE')
        else:
            #adding offset to the VCSS fluxes:
            offset = 1.075
            flux_arr[my_ind] = offset*flux_arr[my_ind]
            eflux_arr[my_ind] = offset*flux_arr[my_ind]
            
    return freq_arr, flux_arr, eflux_arr, labels
 
def mod_data_flux_density_scale(fluxes, efluxes, labels ):
    '''
    Recommendations by the referee:
        RACS and NVSS: add a measurement uncertainity of 5% in quadrature. 
        Add LOTSS-DR1
        WENSS reduce flux by 10%
        TGSS - increase the errors  by 10%
    '''
    for ii, label in enumerate(labels):
        if label not in [ 'VLASS', 'WENSS']:
            efluxes[ii] = np.sqrt(efluxes[ii]**2 + (0.05*fluxes[ii])**2)
        if label in ['WENSS']:
            fluxes[ii] = 0.81*fluxes[ii]
            efluxes[ii] = 0.81*efluxes[ii] 
            efluxes[ii] = np.sqrt(efluxes[ii]**2 + (0.05*fluxes[ii])**2)
        
    return fluxes, efluxes, labels


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
    
'''vla_ax.add_index('Source_name')
ind = vla_ax.loc_indices['J0643-27']
vla_ax.remove_row(ind)   
vla_ax.remove_indices('Source_name')'''

vla_ax_grp = vla_ax.group_by('WISEname')
vla_bx_grp = vla_bx.group_by('WISEname')
    
sp_class = loadt.get_sp_info()
sp_class_sources = list(sp_class[sp_class['Obs_Flag']!='NO']['Source_name'])
sp_class.add_index('Source_name')

spdf = pd.read_excel('New_Analysis/Spectral_shape_classf_master.xlsx')
sp_class = Table.from_pandas(spdf)
sp_class.add_index('Source_name')

#additional VLITE data

vlite_sources = {'J1238+52':17,  'J0823-06':16.2, 'J0642-27':20 }

# fitTab = ascii.read('GuessPar_Radiofits.csv', format='basic', delimiter=',')
# fitTab.add_index('Name')

fit_results_ext_tab = Table(names = ('source', 'region',  'PL_s0', 'PL_al', 'CPL_s0', 
                                'CPL_al','CPL_nup',  'EFFA_s0', 'EFFA_al', 'EFFA_nup',
                                 'IFFA_s0', 'IFFA_al', 'IFFA_nup','SSA_s0', 
                                 'SSA_al', 'SSA_nup', 'GCV_s0', 'GCV_al1', 'GCV_al2', 
                                 'GCV_nup', 'perr_pl_s0', 'perr_pl_al',
                                 'perr_cpl_s0','perr_cpl_al', 'perr_cpl_nup', 
                                 'perr_effa_s0','perr_effa_al', 'perr_effa_nup',
                                 'perr_iffa_s0','perr_iffa_al', 'perr_iffa_nup',
                                 'perr_ssa_s0','perr_ssa_al', 'perr_ssa_nup',
                                 'perr_gcv_s0', 'perr_gcv_al1', 'perr_gcv_al2', 
                                 'perr_gcv_nup',
                                 'rchisq_pl', 'rchisq_cpl','rchisq_effa',
                                 'rchisq_iffa','rchisq_ssa', 'rchisq_gcv'
                                 ), 
                        dtype = ('S10', 'S10', 'f8', 'f8', 'f8', 'f8', 'f8',
                                 'f8', 'f8', 'f8','f8', 'f8', 'f8','f8', 'f8', 
                                 'f8','f8','f8','f8','f8','f8','f8','f8','f8',
                                 'f8','f8','f8','f8','f8','f8','f8','f8','f8',
                                 'f8','f8', 'f8', 'f8','f8', 'f8',
                                 'f8', 'f8','f8', 'f8', 'f8'))
    
cpl_fits = dict.fromkeys(wise_phot['Source_name'])
for key, value in cpl_fits.items():
    cpl_fits[key] = dict.fromkeys(['s0_cpl', 'al_cpl', 'q_cpl', 'speak', 'nupeak' ])
    

mcsim_flag = True
wise_phot.sort('WISEname')

useInBand = True

myfluxdict = dict.fromkeys(list(wise_phot['Source_name']))

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
    morph = sp_row['Xband_morph']
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
    OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,
                                           jvla_BX)
    
    freq_arr, flux_arr, eflux_arr, labels = mcsim_fit.modify_data(sname,freq_arr, flux_arr,
                                eflux_arr, list(labels),sp_row, scat  )
    
    #get VCSS data
    #freq_arr, flux_arr, eflux_arr, labels = get_VCSS_data(sname, wise_row,
    #                         sp_row, labels, freq_arr, flux_arr, eflux_arr)
    freq_arr, flux_arr, eflux_arr, labels = VCSS_data_cleanup(sname, wise_row,
                            sp_row, labels, freq_arr, flux_arr, eflux_arr)
    
    flux_arr, eflux_arr, labels = mod_data_flux_density_scale(flux_arr, eflux_arr, 
                                                              labels)
    
    if 'VLASS' in labels and sp_row['VLASS_Limit']!='U':
        ind = labels.index('VLASS')
        vflux = flux_arr[ind]
        eflux = eflux_arr[ind]
        if eflux<0.2*vflux:
            eflux_arr[ind] = np.sqrt(eflux**2 + (vflux/5)**2)
            
    if useInBand:
        if 'BX' in labels and np.any(jvla_BX['snr']>50):
            freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, 
                                            flux_arr, eflux_arr, alpha_BX)
        if 'AX' in labels and np.any(jvla_AX['snr']>50):
            freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, 
                                            flux_arr, eflux_arr, alpha_AX)

    #create dict for online table:
    mydict = {}
    for snu, esnu, cat in zip(flux_arr, eflux_arr, labels):
        mydict['F'+cat] = np.round(snu,2)
        mydict['E'+cat] = np.round(esnu,2)
        
    if 'BX' in labels and np.any(jvla_BX['snr']>50) and alpha_BX!=[]:
        if alpha_BX[1]>-10:
            mydict['SI_BX'] = alpha_BX[1]
            mydict['ESI_BX'] = alpha_BX[2]
    if 'AX' in labels and np.any(jvla_AX['snr']>50) and alpha_AX!=[]:
        if alpha_AX[1]>-10:
            mydict['SI_AX'] = alpha_AX[1]
            mydict['ESI_AX'] = alpha_AX[2]  
        
    myfluxdict[sname] = mydict
    # Guess parameters:
    #guess_pars = fitTab.loc[wisen]
    models = ['PL', 'CPL', 'EFFA', 'IFFA', 'SSA', 'GCV']
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
 
    
with open('cpl_fits.json','w') as fp:
    json.dump(cpl_fits, fp)    

fit_results_ext_tab.write('Fit_results_new_VCSS_RACS_inband_fscale.csv', overwrite= True)
#fit_results_ext_tab.write('Fit_results_new_VCSS_RACS.csv', overwrite= True)

mydf = pd.DataFrame.from_dict(myfluxdict,   orient='index')
mydf.to_csv('Tab4_Flux_measurements_spectra_fscale.csv')

'''Following not needed if repeating the SED fitting'''
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

####Adding VCSS limits:
# Change: Feb 12 2021

wise_phot.add_column(Column(data = [np.nan]*len(wise_phot), name = 'VLITE_Limit',
                            dtype='f8'))
vcss_cutout_tab = Table.read('Astroquery/VCSS/Pallavi_vcss_mosaics_filenames_5arcsec_all.fits')
vcss_det_tab = Table.read('../Radio-SED-Fitting/Datasets/pallavi_vcss_5arcsec_all.fits')

vcss_det_tab.add_index('WISEname')
vcss_cutout_tab.add_index('WISEname')


for ii, row in enumerate(wise_phot):
    wisen = row['WISEname']
    vcss_det = vcss_det_tab.loc[wisen]
    vcss_msc = vcss_cutout_tab.loc[wisen]
    if vcss_det['index_1']>0:
        continue
    else:
        if vcss_msc['noise']>0:
            wise_phot[ii]['VLITE_Limit'] = vcss_msc['noise']
    
    
wise_phot.write('WISE_phot_comb_manual_fluxes.csv') 




final_unq_jvla = Table.read('New_Analysis/final_unq_JVLASP_vf_morph_regmod.csv')
final_unq_jvla.add_index('Source_name')

fp = open('Flux_info_dict.json', 'r')
flux_dict = json.load(fp)
fp.close()

# Create a dictionary that will specify the list of catalogs used in the 
# calculation of alpha_high and alpha_low


def calc_nu_p(alpha,q ):
    return np.exp(-alpha/(2*q))

survey_info = dict.fromkeys(['NVSS', 'FIRST', 'SUMSS', 'TGSS', 'VLSSr', 'WENSS',
                          'WISH', 'GB6', 'VLITE', 'VLASS', 'TEXAS',  'GLEAM', 'LOTSS', 
                          'AX', 'BX'])
freq_survey = [1.4, 1.4, 0.843, 0.150, 0.074, 0.325, 0.352, 4.85, 0.338, 3.0, 0.320, 
                           0.200, 0.150, 10, 10]

for key, freq in zip(survey_info.keys(), freq_survey):
    fdict = {'Freq':freq}
    survey_info[key] = fdict

'''
Calculation of alhigh and allow using different surveys
'''
al_calcs  = Table(names=('Source_name', 'al_l_TGSS_NVSS', 'e_al_l_TGSS_NVSS',
                'l_al_l_TGSS_NVSS', 'al_h_NVSS_VLAX', 'e_al_h_NVSS_VLAX',
                'l_al_h_NVSS_VLAX', 'al_l_WENSS_NVSS', 'e_al_l_WENSS_NVSS',
                'l_al_l_WENSS_NVSS', 'al_h_VLASS_VLAX', 'e_al_h_VLASS_VLAX',
                'l_al_h_VLASS_VLAX'), dtype=('S10', 'f8','f8','bool', 'f8', 'f8',
                                   'bool', 'f8', 'f8', 'bool', 'f8', 'f8', 'bool'))

for source in final_unq_jvla['Source_name']:
    row= [str(source)]
    sub = [np.nan, np.nan, False]
    row.extend(sub*4)
    al_calcs.add_row(row)
    
for ii, row in enumerate(al_calcs):
    source = row['Source_name']
    sp_row = sp_class.loc[source]
    pardict = flux_dict[source]
    labels = list(pardict.keys())
    freq_arr = np.array([pardict[k][0] for k in labels])
    flux_arr = np.array([pardict[k][1] for k in labels])
    eflux_arr = np.array([pardict[k][2] for k in labels])

    hcat = 'BX'
    if 'BX' not in labels:
        hcat = 'AX'
    lcat = 'NVSS'
    nu =  [pardict[hcat][0], pardict[lcat][0]]
    s  =  [pardict[hcat][1], pardict[lcat][1]]
    e_s = [pardict[hcat][2], pardict[lcat][2]]
        
    sp, e_sp = dquants.two_pt_alpha(nu, s, e_s)
    al_calcs[ii]['al_h_NVSS_VLAX', 'e_al_h_NVSS_VLAX'] = sp, e_sp
    
    
    
    if 'VLASS' in labels:
        lcat = 'VLASS'
        nu =  [pardict[hcat][0], pardict[lcat][0]]
        s  =  [pardict[hcat][1], pardict[lcat][1]]
        e_s = [pardict[hcat][2], pardict[lcat][2]]
        if ~np.isnan(s).any():  
            sp, e_sp = dquants.two_pt_alpha(nu, s, e_s)
            al_calcs[ii]['al_h_VLASS_VLAX', 'e_al_h_VLASS_VLAX'] = sp, e_sp
            lcol = lcat+'_Limit'
            if lcol in sp_row.colnames:
                if sp_row[lcol] == 'U':
                    al_calcs[ii]['l_al_h_VLASS_VLAX'] = True
    
    hcat = 'NVSS'
    lcat = 'TGSS'
    
    if 'TGSS' in labels:
        nu =  [pardict[hcat][0], pardict[lcat][0]]
        s  =  [pardict[hcat][1], pardict[lcat][1]]
        e_s = [pardict[hcat][2], pardict[lcat][2]]
        if ~np.isnan(s).any():  
            sp, e_sp = dquants.two_pt_alpha(nu, s, e_s)
            al_calcs[ii]['al_l_TGSS_NVSS', 'e_al_l_TGSS_NVSS'] = sp, e_sp
            lcol = lcat+'_Limit'
            if lcol in sp_row.colnames:
                if sp_row[lcol] == 'U':
                    al_calcs[ii]['l_al_l_TGSS_NVSS'] = True

    hcat = 'NVSS'
    lcat = 'WENSS'
    
    if 'WENSS' in labels:
        nu =  [pardict[hcat][0], pardict[lcat][0]]
        s  =  [pardict[hcat][1], pardict[lcat][1]]
        e_s = [pardict[hcat][2], pardict[lcat][2]]
        if ~np.isnan(s).any():  
            sp, e_sp = dquants.two_pt_alpha(nu, s, e_s)
            al_calcs[ii]['al_l_WENSS_NVSS', 'e_al_l_WENSS_NVSS'] = sp, e_sp
            lcol = lcat+'_Limit'
            if lcol in sp_row.colnames:
                if sp_row[lcol] == 'U':
                    al_calcs[ii]['l_al_l_WENSS_NVSS'] = True

     
    
    
al_calcs = join(al_calcs, sp_class, keys= 'Source_name', join_type='left')        
al_calcs.add_index('Source_name')      


sed_fit_params= Table(names =('Source_name', 'SED_Class', 'al_10_1.4', 'e_al_10_1.4', 'l_al_10_1.4', 
                             'al_1.4_.15', 'e_al_1.4_.15', 'l_al_1.4_.15', 'nu_p', 'e_nu_p', 'l_nu_p',
                              's_p', 'e_s_p', 'q', 'e_q', 'al_highc', 'e_al_highc', 'l_al_highc', 
                             'al_lowc', 'e_al_lowc', 'l_al_lowc','quality', 'ccode', 's0_highc', 
                             'e_s0_highc', 's0_lowc', 'e_s0_lowc', 'e_nu_p_eh', 'e_nu_p_el',
                             'e_s_p_eh', 'e_s_p_el' ), 
                         dtype = ('S8', 'S3', 'f8', 'f8', 'bool', 'f8', 'f8', 'bool','f8', 'f8', 'bool',
                                 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'bool','f8', 'f8', 'bool',
                                 'S2', 'S1', 'f8', 'f8', 'f8', 'f8','f8', 'f8', 'f8', 'f8'))


    
fit_results_ext_tab.add_index('source')

for row in final_unq_jvla:
    source = row['Source_name']
    al_row = al_calcs.loc[source]
    spshape =  al_row['Final_SP_Class']   
    fitRes = fit_results_ext_tab.loc[source]
    par_cpl = fitRes['CPL_s0', 'CPL_al', 'CPL_nup']
    perr_cpl = fitRes['perr_cpl_s0','perr_cpl_al','perr_cpl_nup']
    nu_p = calc_nu_p(par_cpl[1],par_cpl[2])
    s_p = rmfit.CPL(nu_p, *par_cpl)
    
    myrow = [source, al_row['Final_SP_Class'], al_row['al_h_NVSS_VLAX'], al_row['e_al_h_NVSS_VLAX'], 
            al_row['l_al_h_NVSS_VLAX'], al_row['al_l_TGSS_NVSS'], al_row['e_al_l_TGSS_NVSS'], 
            al_row['l_al_l_TGSS_NVSS'], np.nan, np.nan, False, np.nan, np.nan, par_cpl[2], perr_cpl[2],
            np.nan, np.nan, False, np.nan, np.nan, False, al_row['Sp_qual'], al_row['Clean_Code'], 
            np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    
    if spshape in ['CV', 'CV:', 'PK', 'PK:', 'I']:
        
        myrow[8] = nu_p
        myrow[11] = s_p

        e_nup_fit = emcsim.errors_mcsim(calc_nu_p, 2, 2, [par_cpl[1], par_cpl[2]], [perr_cpl[1], perr_cpl[2]],
                                                            10000, True, False)
        myrow[9] = e_nup_fit['Std dev']
        nup_err = e_nup_fit['Std dev']
        if nup_err>1000:
            err1 = np.abs(e_nup_fit['e1High']-e_nup_fit['Central Value'])
            err2 = np.abs(e_nup_fit['e1Low']-e_nup_fit['Central Value'])
            nup_err = np.mean([err1, err2])
            myrow[27] = err1
            myrow[28] = err2
            
        e_sp_fit =  emcsim.errors_mcsim(rmfit.CPL, 4, 4, [nu_p, par_cpl[0], par_cpl[1],par_cpl[2]],
                            [nup_err,perr_cpl[0], perr_cpl[1],perr_cpl[2]], 10000, True, False)
        
        myrow[12] = e_sp_fit['Std dev']
        if e_sp_fit['Std dev']>1000:
            err1 = np.abs(e_sp_fit['e1High']-e_sp_fit['Central Value'])
            err2 = np.abs(e_sp_fit['e1Low']-e_sp_fit['Central Value'])
            myrow[29] = err1
            myrow[30] = err2
            
        print(nu_p)    
    sed_fit_params.add_row(myrow)
        
sed_fit_params.add_index('Source_name')



aldict = dict.fromkeys(list(final_unq_jvla['Source_name']))
for k, v in aldict.items():
    aldict[k] = {'al_high':None, 'al_low':None}
for row in final_unq_jvla:
    source = row['Source_name']
    al_row = al_calcs.loc[source]
    sp = al_row['Final_SP_Class']
    aldict[source]['SpShape'] = sp
    
for row in final_unq_jvla:
    source = row['Source_name']
    kk = sed_fit_params.loc_indices[source] 
    wise_row = wise_phot.loc['Source_name', source]
    wisen = wise_row['WISEname']
    scat = wise_row
    glmcat = atgcat.loc[wisen]
    al_row = al_calcs.loc[source]
    spshape = al_row['Final_SP_Class']
    sp_row = sp_class.loc[source]
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
    OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,
                                           jvla_BX)
    freq_arr, flux_arr, eflux_arr, labels = mcsim_fit.modify_data(source,freq_arr, flux_arr,
                                eflux_arr, labels,sp_row, scat  )
    
    if 'VLASS' in labels and sp_row['VLASS_Limit']!='U':
        ind = labels.index('VLASS')
        vflux = flux_arr[ind]
        eflux = eflux_arr[ind]
        if eflux<0.2*vflux:
            eflux_arr[ind] = np.sqrt(eflux**2 + (vflux/5)**2)        

    hdict = []
    ldict = []
    
    if 'BX' in labels and np.any(jvla_BX['snr']>50):
        freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, 
                                        flux_arr, eflux_arr, alpha_BX)
        hdict.append('IB')
    if 'AX' in labels and np.any(jvla_AX['snr']>50):
        freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, 
                                        flux_arr, eflux_arr, alpha_AX)
        hdict.append('IB')

    msk = [freq_arr>0]
    nu_arr = freq_arr[tuple(msk)]
    s_arr= flux_arr[tuple(msk)]
    es_arr = eflux_arr[tuple(msk)]

    
    #calculation of alhigh and alow:
    #Strategy: Identify peak and calculate points before and after peak location
    min_cat = labels[np.argmin(nu_arr)]
    min_nu = survey_info[min_cat]['Freq']
    max_cat = 'BX'
    max_nu = 10
    if 'AX' in labels:
        max_cat = 'AX'
    nup = sed_fit_params[kk]['nu_p']
    e_nup = sed_fit_params[kk]['e_nu_p']
    if nup<10 and nup>min_nu:
        hind = np.nonzero(nu_arr>nu_p)[0]
        lind = np.nonzero(nu_arr<nu_p)[0]        
        for ind in hind:
            hdict.append(labels[ind])
        
        for ind in lind:
            ldict.append(labels[ind])         
        aldict[source]['al_high'] = hdict
        aldict[source]['al_low'] = ldict
        
    if nu_p>0 and nup<min_nu:
        hind = np.nonzero(nu_arr>1.3)[0]
        for ind in hind:
            hdict.append(labels[ind])
        aldict[source]['al_high'] = hdict

    if spshape=='I':
        for ind in range(len(labels)):
            ldict.append(labels[ind])
        aldict[source]['al_low'] = ldict

    if spshape=='PL' or spshape=='PL:':
        for ind in range(len(labels)):
            hdict.append(labels[ind])
        aldict[source]['al_high'] = hdict
    
#Writing the dictionary:
with open('Alpha_calc_info_v2.json', 'w') as fp:
    json.dump(aldict, fp, indent=4)