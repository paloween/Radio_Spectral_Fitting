#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 22:48:59 2018

@author: pallavipatil
"""

# Estimating Space densities for WISE color selection. 

import glob
import math
import os
import numpy as np
import astropy.io.fits as fits
from astropy.visualization import ZScaleInterval
import aplpy
import montage_wrapper as montage
import matplotlib
from astropy.table import Table, join,unique, Column,hstack,vstack
from matplotlib import pyplot as plt
from astropy.io import fits,ascii
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
import astropy.units as u
#from file_prep_radio_fitting import data_prep
import aplpy as apl
from astropy.visualization import ZScaleInterval
import scipy
import pandas as pd
#import Radio_Models_func as rm
import astropy.constants as const
from astropy.wcs import WCS
from astropy.modeling import models, fitting
import matplotlib.patches as mpatches
import matplotlib.cm as cm 
from astropy.coordinates import Angle
from astropy.visualization.wcsaxes import SphericalCircle
import matplotlib.hatch as mhatch
from astropy.cosmology import Planck15 as cml
from scipy.interpolate import interp1d

# WISE AllSky Catalog 
wise_allsky = Table.read('/Users/pallavipatil/Desktop/VLA/WISE_AllSky.fits')
wise_allsky.rename_column('wisename_01', 'WISEname')
wise_phot = Table.read('../WISE-photometry.csv')
wise_sname = []
for wisen in wise_phot['WISEname']:
    ra= wisen[:4]
    if len(wisen.split('-'))>1:
        sign = '-'
    else:
        sign= '+'
    dec = wisen.split(sign)[1][:2]
    wise_sname.append('J'+ra+sign+dec)
wise_phot.add_column(Column(data=wise_sname, name = 'Source_name'))
wise_allsky.add_column(Column(data=wise_sname, name = 'Source_name'))
wise_phot = join(wise_phot, wise_allsky, keys='Source_name' )
final_unq_jvla = Table.read('final_unique_jvla_tab.csv', format = 'csv',)
final_unq_jvla.add_index('Source_name')
wise_phot.add_index('Source_name')
# Groups based on morphology 


# WISE bands
WISE_l = np.array([3.35, 4.60, 11.56, 22.08])
WISE_bw = np.array([0.66, 1.04, 5.51, 4.10])
WISE_nu = const.c.value*10**6/np.array(WISE_l)
WISE_name = np.array(['W1', 'W2', 'W3', 'W4'])
WISE_tab = Table(data = (WISE_name, WISE_l, WISE_bw), names = ('Band','Wlam', 'BW'))
lam_obs = [3.35, 4.60, 11.56, 22.08]
lam_obs.append(214130)
lam_obs = np.array(lam_obs)

# WISE Sample selection
def WISE_color_crit(W1, W2, W3):
    return W1-W2 + 1.25*(W2-W3) - 7

#WISE SNR selection
def WISE_snr_selection(W3snr, W4snr):
    if W3snr > 7 or W4snr > 7:
        return True
    else:
        return False
    
#WISE W4 flux selection
def WISE_W4_flux_crit(W4):
    if W4>7: # In mJy
        return True
    else:
        return False
    
    
    
def check_for_linefeat(z):
    # PAH & silicate features:
    PAH_feat = [6.2, 7.7, 8.6, 11.3]
    silicate_feat = [9.7, 18]
    all_feat = PAH_feat
    all_feat.extend(silicate_feat)
    WISE_l = np.array([3.35, 4.60, 11.56, 22.08])
    WISE_bw = np.array([0.66, 1.04, 5.51, 4.10])
    WISE_rf_l = WISE_l/(1+z)
    WISE_rf_bw = WISE_bw/(1+z)
    WISE_l_low = WISE_rf_l-WISE_rf_bw/2 
    WISE_l_up = WISE_rf_l+WISE_rf_bw/2 
    wise_lfeat = [False, False, False, False]
    for ii in range(0,4):
        l_range = [WISE_l_low[ii], WISE_l_up[ii]]
        for lft in all_feat:
            if lft>l_range[0] and lft <l_range[1]:
                wise_lfeat[ii] = True
    return wise_lfeat
        


def calculate_wise_mags(wise_flux):
    F_nu0  = np.array([306.682, 170.663, 29.045, 8.284])
    fc_al2 = np.array([1.0, 1.0, 1.0, 1.0])
    fc = fc_al2
    mags = -2.5*np.log10((wise_flux*fc*0.001)/F_nu0)
    return mags
    
def radio_flux_calc(fnvss,zz, alpha, zT):
    Dl_T = cml.luminosity_distance(zT).value
    Dl_p = cml.luminosity_distance(zz).value
    snew = (Dl_T/Dl_p)**2*((1+zz)/(1+zT))**(1+alpha)*fnvss
    #nu_rf_zT = 1.4*(1+zT)
    #nu_rf_z = 1.4*(1+z)
    #snew = fnvss*(nu_rf_z/nu_rf_zT)**(alpha)
    return snew

def sample_selection(wise_flux, W3snr, W4snr, wise_mags,z, fnvss_z, q,fnvss):
    reason = []
    result_flag = True
    W4_flux = WISE_W4_flux_crit(wise_flux[3])
    '''if W4_flux==False: 
        result_flag = False
        reason.append(['W4_flux'])'''
    snr_crit = WISE_snr_selection(W3snr, W4snr)
    if snr_crit == False:
        result_flag = False
        reason.append(['W3-W4_SNR'])
    '''if wise_flux[2] < W3lim:
        result_flag = False
    if wise_flux[3] < W4lim:
        result_flag = False'''
    wise_color_crit = WISE_color_crit(wise_mags[0], wise_mags[1], wise_mags[2])
    if wise_color_crit<0:
        result_flag = False 
        reason.append(['WISE_colors'])
    if fnvss > 7:
        if fnvss_z < 7:
            result_flag = False
            reason.append(['NVSS_flux'])
    if q>0:
        result_flag = False
        reason.append(['q-RQ'])       
    return result_flag, reason

def fit_poly(x, y, deg, xE):
    f1 = np.polyfit(x, y, deg)
    func = np.poly1d(f1)
    return func(xE)

def WISE_flux_intpl(z, wise_flux, WISE_l, WISE_nu, zT):        
    #Rest frame Lam at True z:
    wise_nurf_T = WISE_nu*(1+zT)
    wise_nurf = WISE_nu*(1+zz)
    wise_lrf = WISE_l/(1+zz)
    
    #Luminosity distance
    Dl_T = cml.luminosity_distance(zT).value
    Dl_p = cml.luminosity_distance(zz).value
    #calculating spectral luminosities
    L_nu_rf = 4*np.pi*Dl_T*Dl_T*np.array(list(wise_flux))/(1+zT)

    #working in log-log space
    log_nu = np.log10(wise_nurf_T)
    log_Lnu = np.log10(L_nu_rf)
    log_nu_p = np.log10(wise_nurf)
    fit_Lnu = np.polyfit(log_nu, log_Lnu, 2)
    eval_Lnu = np.poly1d(fit_Lnu)
    
    #evaluating Luminosities
    Lnu_p = 10**(eval_Lnu(log_nu_p))
    #evaluating fluxes
    wise_flux_p = Lnu_p *(1+zz)/(4*np.pi*Dl_p*Dl_p)    
    return wise_flux_p, wise_nurf, wise_lrf
    '''
    wise_alpha = scipy.stats.linregress(np.log10(WISE_nu), np.log10(list(wise_flux))).slope
    f_W3 = wise_flux[2]
    nu_W3_rf = (WISE_nu*(1+zT))[2]
    #assuming f = nu**(wise_alpha)
    #f/f(W3) = (nu/nu(w3))**wise_alpha
    wise_fz = f_W3*(wise_nurf/nu_W3_rf)**wise_alpha
    return wise_fz, wise_nurf, wise_lrf'''



vmax_arr = []
vmax_L = []
vmax_shell = []
vmax_z = []
v_z = []
zmax_arr = []

for name in final_unq_jvla['Source_name']:
    z = wise_phot.loc[name]['redshift']
    if z>0:
        vmax_z.append(z)
plt.ioff()

#Note: calculate the vmax for sources without z
#set the zflag to flase to skip this calculation
zflag = False
if zflag:
    for name in final_unq_jvla['Source_name']:
        z = wise_phot.loc[name]['redshift']
        if z>0:
            vmax_z.append(z)
        else:
            rdmz = np.random.uniform(0.4, 2.8)
            vmax_z.append(rdmz)

for sind in range(0,len(final_unq_jvla)):
    #Looping over z:
    zinc = 0.01
    zmin = 0.01
    zmax = 6.0
    nzs = int((zmax-zmin)/zinc)
    zarr = np.linspace(0, 6.0, nzs)
    Asky_common = 28443 #Sq deg
    Atot = 4*np.pi*(180/np.pi)**2
    SNR_cut = 7 
    jvla_tab = final_unq_jvla[sind]
    sname = jvla_tab['Source_name']

    row = wise_phot.loc[sname]


    wise_flux_colnames = ['FWISE1','FWISE2','FWISE3','FWISE4']
    wise_mag_colnames = ['w1mpro', 'w2mpro', 'w3mpro', 'w4mpro']
    W1ind = wise_phot.index_column('FWISE1')
    wise_flux = row[wise_flux_colnames]
    zT = row['redshift']
    if zflag and not zT>0: zT = vmax_z[sind]
    fnvss = row['NVSS1.4']
    if zT>0:
        if jvla_tab['alh_BX'] >-99:
            alpha = jvla_tab['alh_BX']
        elif jvla_tab['alh_AX'] >-99:
            alpha = jvla_tab['alh_AX']
        else:
            alpha = -0.9    


        W3rms = wise_flux['FWISE3']/row['w3snr']
        W4rms = wise_flux['FWISE4']/row['w4snr']
        wise_mags = row[wise_mag_colnames]
        W3lim = SNR_cut*W3rms
        W4lim = SNR_cut*W4rms
        W1_W2 = wise_mags[0]- wise_mags[1]
        W2_W3 = wise_mags[1]- wise_mags[2]
        q = np.log10(wise_flux['FWISE4']/fnvss)
        
        zspace = np.arange(zmin, zmax, zinc)
        info_row = [sname, zT,True]
        info_row.extend(WISE_l/(1+zT))
        info_row.extend(list(wise_flux))
        info_row.extend([row['w3snr'],row['w4snr']])
        info_row.extend(wise_mags)
        info_row.extend([W1_W2, W2_W3])
        info_row.extend([False, False, False, False])
        info_row.extend([fnvss,q,''])
        info_names = ['Source_name', 'z_inc', 'Sample_Selection', 'W1_rf', 'W2_rf','W3_rf',
                      'W4_rf', 'FW1_rf','FW2_rf','FW3_rf','FW4_rf',
                      'W3snr', 'W4snr', 'mW1_rf','mW2_rf','mW3_rf','mW4_rf',
                      'W1-W2', 'W2-W3','W1_lft', 'W2_lft','W3_lft','W4_lft','FNVSS', 'q', 'Reason']
        wise_vmax_calc = Table( names = info_names, dtype = (['S10','f8','S1']+[ 'f8']*16+['S1']*4+['f8','f8']+['S30']))
        wise_vmax_calc.add_row(info_row)
        vmax_s = 0 
        
        wise_lrf_T = WISE_l/(1+zT)
        wise_tmp_lam = np.logspace(-3,3, 100)
        wise_tmp = np.poly1d(np.polyfit(np.log10(wise_lrf_T), np.log10(list(wise_flux)), deg = 2))
        
        
        for zz in zspace:
            wise_nurf = np.array(WISE_nu)*(1+zz)
            wise_lrf = np.array(WISE_l)/(1+zz)
            #wise_fz = 10**(wise_tmp(np.log10(wise_lrf)))
            
            wise_fz, wise_nurf, wise_lrf = WISE_flux_intpl(zz, wise_flux, WISE_l, WISE_nu, zT)
            # Check for sample selections
            wise_magsz = calculate_wise_mags(wise_fz)
            w3snr = wise_fz[2]/W3rms
            w4snr = wise_fz[3]/W4rms
            fnvss_z = radio_flux_calc(fnvss, zz, alpha, zT)
            q = np.log10(wise_fz[3]/fnvss_z)
            wise_sample_z,reason = sample_selection(wise_fz,w3snr, w4snr, wise_magsz, zz, fnvss_z, q,fnvss )
            wise_linft = check_for_linefeat(zz)
            w1_w2 = wise_magsz[0] - wise_magsz[1]
            w2_w3 = wise_magsz[1] - wise_magsz[2]

            #print(wise_magsz, w1_w2, w2_w3)
            #print(wise_fz, wise_magsz, wise_sample_z, wise_linft)
            info_row = [sname, zz,wise_sample_z ]
            info_row.extend(wise_lrf)
            info_row.extend(wise_fz)
            info_row.extend([w3snr, w4snr])
            info_row.extend(wise_magsz)
            info_row.extend([w1_w2, w2_w3])
            info_row.extend(wise_linft)
            info_row.extend([fnvss_z, q, str(reason)])
            wise_vmax_calc.add_row(info_row)
            if wise_sample_z:
                vol_z = cml.comoving_volume(zz+zinc) - cml.comoving_volume(zz)
                vmax_s = vmax_s + vol_z
        
        vmax = 0
        zmax_iarr = np.where(wise_vmax_calc[wise_vmax_calc['z_inc']>zT]['Sample_Selection']=='F')
        if len(zmax_iarr[0]) >0:
                zmax_ind = zmax_iarr[0]
                zmax = wise_vmax_calc[wise_vmax_calc['z_inc']>zT]['z_inc'][zmax_ind][0]
                dlmax = cml.luminosity_distance(zmax).value
                vmax = cml.comoving_volume(zmax).value
        else:
            zmax = 6.0
            dlmax = cml.luminosity_distance(zmax).value
            vmax = cml.comoving_volume(zmax).value
        zmax_arr.append(zmax)    
        print(zmax, zT)
        covol = cml.comoving_volume(zT).value
        v_z.append(covol)
        vmax_arr.append(vmax)
        vmax_shell.append(vmax_s)
        dl = cml.luminosity_distance(zT).to('m').value
        Lnvss = 4*np.pi*dl*dl*fnvss*10**(-3)*10**(-26)/(1+zT)**(1+alpha)
        vmax_L.append(np.log10(Lnvss))
        ascii.write(wise_vmax_calc, './Vmax_calc/'+sname+'_vmax_tab.csv', format = 'csv',
                        delimiter = ',', overwrite = True)
        print(sname, zmax, "{:.2e}".format(vmax),"{:.2e}".format(vmax_s) )
        
        #  PART to create SED plots:
        zmin = 0.01
        zmax = 6.00

        wise_lrf_lz = np.array(lam_obs)/(1+zmin)
        lf_lz = list(wise_vmax_calc['FW1_rf','FW2_rf','FW3_rf','FW4_rf'][1])
        #lf_lz = 10**(wise_tmp(np.log10(wise_lrf_lz[:-1])))
        fnvss_lz = radio_flux_calc(fnvss, zmin, alpha, zT)
        lf_lz = np.append(lf_lz, fnvss_lz)
        wise_lrf_hz = np.array(lam_obs)/(1+zmax)
        #lf_hz = 10**(wise_tmp(np.log10(wise_lrf_hz[:-1])))
        lf_hz = list(wise_vmax_calc['FW1_rf','FW2_rf','FW3_rf','FW4_rf'][-1])
        fnvss_hz = radio_flux_calc(fnvss, zmax, alpha, zT)
        lf_hz = np.append(lf_hz, fnvss_hz)
        wise_tmp = np.poly1d(np.polyfit(np.log10(wise_lrf_T), np.log10(list(wise_flux)), deg = 2))

        flux_true = list(wise_flux)
        flux_true.append(fnvss)
        lam_true = lam_obs
        ftmp = 10**wise_tmp(np.log10(wise_tmp_lam))
        fig, axs = plt.subplots(2,1)
        ax = axs[0]
        ax.loglog(lam_true/(1+zT), flux_true, 'o', label='Observed', color='blue')
        ax.loglog(wise_lrf_lz, lf_lz, '*', label='z=0.01', color='green')
        ax.loglog(wise_lrf_hz, lf_hz, '^', label='z=6', color= 'red')
        ax.loglog(wise_tmp_lam, ftmp, ls = 'dotted', label='Fitted Tmp', color= 'grey', alpha = 0.5)
        ax.set_ylim(10**-7,)
        ax.invert_xaxis()
        ax.set_xlabel('Rest Frame, '+r'$\lambda (\mu)$'+'m' )
        ax.set_ylabel('Flux, mJy' )
        ax.legend()
        
        ax1 = axs[1]
        ax1.loglog(lam_obs, flux_true, 'o', label='zT='+str(zT), color='blue')
        ax1.loglog(wise_tmp_lam*(1+zT), ftmp, ls = 'dotted', alpha = 0.5, color='blue')
        ax1.loglog(lam_obs, lf_lz, '*', label='z=0.01', color='green')
        ax1.loglog(wise_tmp_lam*(1+zmin), ftmp, ls = 'dotted', alpha = 0.5, color='green')
        ax1.loglog(lam_obs, lf_hz, '^', label='z=6.0', color='red')
        ax1.loglog(wise_tmp_lam*(1+zmax), ftmp, ls = 'dotted', alpha = 0.5, color='red')
        ax1.set_ylim(10**-7, )
        ax1.invert_xaxis()
        ax1.set_xlabel('Observed, '+r'$\lambda (\mu)$'+'m' )
        ax1.set_ylabel('Flux, mJy' )
        ax1.legend()
        plt.tight_layout()
        plt.savefig('Vmax_SED/'+sname+'_vmax_SED.png')
        plt.close()

        

        
        
            
#plt.scatter(wise_vmax_calc['W1-W2'], wise_vmax_calc['W2-W3'])

#DIagnostics:
        
files= glob.glob('Vmax_calc/*.csv')
plt.ioff()
#os.mkdir('Vmax_plots')
for fi in files:   
    #tab = Table.read('Vmax_calc/J0010+16_vmax_tab.csv')
    tab = Table.read(fi)[1:]    
    fig, axs = plt.subplots(3, 2, figsize=(7,10))
    ax = axs[0,0]
    
    truth_arr = tab['Sample_Selection']
    for ii, t in enumerate(truth_arr):
        if t == 'F':
            truth_arr[ii] = 0
        else:
            truth_arr[ii] = 1
            
            
    ax.plot(np.array(tab['z_inc'])[1:],np.array(truth_arr)[1:])
    
    ax.set_xlabel('Redshift, '+r'$z$')
    ax.set_ylabel('Sample Selection')
    
    
    ax = axs[1,0]
    ax.semilogy(tab['z_inc'], tab['FNVSS'], '*', markersize = 1)
    ax.axhline(7, ls='dotted', color='k')
    ax.set_xlabel('Redshift, '+r'$z$')
    ax.set_ylabel('NVSS Flux, mJy')
    
    ax = axs[2,0]
    ax.plot(tab['z_inc'], tab['q'], '*', markersize = 1)
    ax.axhline(0, ls='dotted', color='k')
    ax.set_xlabel('Redshift, '+r'$z$')
    ax.set_ylabel('q, mJy')
    
    
    ax = axs[0,1]
    csc = cm.get_cmap('jet')
    sc = ax.scatter(tab['W2-W3'], tab['W1-W2'],c = tab['z_inc'], marker = 'o', cmap=csc)
    x = np.linspace(-1, 7, 50)
    y = 7-x*1.25
    ax.plot(x, y, '--', color = 'k', linewidth = 1)
    cb = plt.colorbar(sc, ax =ax)
    cb.set_label('Redshift, '+r'$z$')
    ax.set_xlabel('W2-W3')
    ax.set_ylabel('W1-W2')
    
    ax = axs[1,1]
    ax.plot(tab['z_inc'],tab['W3snr'],'.', color='r', label='W3')
    ax.plot(tab['z_inc'],tab['W4snr'],'.', color='g', label='W4')
    ax.axhline( 7, ls='dotted', color='grey', alpha= 0.5)
    ax.set_yscale('log')
    ax.set_xlabel('Redshift, '+r'$z$')
    ax.set_ylabel('SNR')
    ax.legend()
    
    ax  = axs[2,1]
    ax.plot(tab['z_inc'], tab['FW4_rf'], '+', color='orange')
    ax.set_xlabel('Redshift, '+r'$z$')
    ax.set_ylabel('W4, mJy')
    ax.axhline(7, ls='dotted', color='blue')
    ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig('Vmax_plots/'+fi[10:-4]+'.png')
    plt.close()
    
    
    # Other plot diagnostics:
'''   
plt.ion()
#os.mkdir('Vmax_plots')
for fi in files:   
    #tab = Table.read('Vmax_calc/J0010+16_vmax_tab.csv')
    tab = Table.read(fi)    
    fig, axs = plt.subplots(3, 2, figsize=(7,10))
    ax = axs[0,0]
    
    truth_arr = tab['Sample_Selection']
    for ii, t in enumerate(truth_arr):
        if t == 'F':
            truth_arr[ii] = 0
        else:
            truth_arr[ii] = 1
            
    
for r in tab[:10]:
     x = list(r['W1_rf','W2_rf','W3_rf','W4_rf'])
     y = list(r['FW1_rf','FW2_rf','FW3_rf','FW4_rf'])
     plt.plot(x, y, 'o',  ms = 2)


'''

# calculation of phi

L_bins = np.arange(24, 29, 0.5)
phi_calc = [0]*len(L_bins)
L_count = [0]*len(L_bins)

phi_arr = [0]*len(L_bins)
phi_sig = [0]*len(L_bins)
vmax_cols = [[]]*len(L_bins)

for vmax, L in zip(vmax_arr, vmax_L):
    ld = np.arange(24, 29, 0.5)
    if L>20 and L<30:
        #ind = np.digitize(L, ld)
        #print(ld, l)
        ld = ld-L
        l = np.argmin(abs(ld))
        lbin_ind = int(l)
        #print(ind, lbin_ind, L)        
        ind = lbin_ind
        phi = (1/0.5)*(1.0/vmax)*Atot/Asky_common
        phi_calc[ind] = phi_calc[ind]+1.0/vmax
        L_count[ind] = L_count[ind]+1
        phi_arr[ind] = phi_arr[ind] + phi
        y = list(vmax_cols[ind])
        y.append(1.0/vmax)
        vmax_cols[ind] = y
        phi_sig[ind]  = phi_sig[ind]+(1.0/vmax)**2
       
        
        '''
        print( L_bins[ind], L)
        phi = (1/0.4)*(1.0/vmax)*Atot/Asky_common
        phi_calc[ind] = phi_calc[ind]+1.0/vmax
        L_count[ind] = L_count[ind]+1
        phi_arr[ind] = phi_arr[ind] + phi
        vmax_cols[ind].append(1.0/vmax)




L_bins
for vmax, L in zip(vmax_arr, vmax_L):
    ld = np.array([23, 24, 25, 26, 27, 28, 29])
    if L>20 and L<30:
        ld = ld-L
        l = np.argmin(abs(ld))
        #print(ld, l)
        lbin_ind = int(l)
        print(lbin_ind, phi_calc[lbin_ind], L_bins[lbin_ind], 1/vmax)
        phi_calc[lbin_ind] = phi_calc[lbin_ind]+1.0/vmax
        phi_sig[lbin_ind]  = phi_sig[lbin_ind]+(1.0/vmax)**2
        vmax_cols[lbin_ind].append(1.0/vmax)
        L_count[lbin_ind] = L_count[lbin_ind]+1
        
        
 '''    

phi_arr = [0]*len(L_bins)
phi_err = [0]*len(L_bins)
for ii,vmax_row in enumerate(vmax_cols):
    nfact = (167/71)
    com_fact = nfact*(1/0.5)*Atot/Asky_common
    phi = np.sum(vmax_row)*com_fact
    phie = np.sqrt(np.sum(np.array(vmax_row)**2))*com_fact
    phi_arr[ii] = phi
    phi_err[ii] = phie
    
        
#fit Broken powerlaw:

def fun_brokenPL(x, amp, x_brk, alpha1, alpha2):
    y = []
    for l in x:
        if l<x_brk:
            y.append(amp*(l/x_brk)**alpha1)
        else:
            y.append(amp*(l/x_brk)**alpha2)
    return np.array(y)
    
    

x = np.array(10**(L_bins))
y = np.array(phi_arr)
yerr = np.array(phi_err)
input_par  = [10**-10, 10**25, -0.4, 0.2]
fit_mod,cov1 = scipy.optimize.curve_fit(fun_brokenPL,  
                                x, y, input_par, yerr)   

    
# from pracy et al 
        
logP1_4 = np.arange(21, 26.5, 0.4)
logphi_all = [-3.58, -3.28, -3.15, -3.55, -4.03, -4.51, -4.95, -5.26, -5.46, -5.65, -6.03, -6.57, -6.83, -7.07]
logphi_agn = [np.nan, np.nan, -3.97, -4.1, -4.56, -4.78, -5.02, -5.28, -5.47, -5.65, -6.03, -6.57, -6.83, -7.07]
logphi_lerg = [np.nan, np.nan, -4.00, -4.19, -4.72, -4.85, -5.09, -5.34, -5.52, -5.71, -6.07, -6.78, -7.26, -7.42]
logphi_herg = [np.nan, np.nan, -4.75, -5.06, -5.07, -5.58, -5.84, -6.15, -6.43, -6.52, -7.04, -6.98, -7.03, -7.33]

# from best 2014
logP_B = np.arange(23.5, 29, 0.5)
best_agn = [-4.67, -5.21, -5.25, -5.41, -6.02, -6.28, -6.73, -7.65, -8.44, -9.37, -10.08]
best_radm = np.array([np.nan, -5.78, -6.34, -6.07, -6.73, -6.46, -6.93, -7.86, -8.52, -9.40, -10.08])
best_radm_err = np.array([[np.nan,np.nan], [0.53, 0.23], [0.70, 0.26], [0.24, 0.15], [0.55, 0.24], 
                 [0.18, 0.12],[0.15, 0.11], [0.10,0.08],[0.08,0.07],[0.13,0.10],[0.30,0.18]])
best_jetm = np.array([-4.79, -5.53, -5.36, -5.74, -6.48, -7.01, -7.39, -8.54,np.nan,np.nan,np.nan ])
best_jetm_err = np.array([[0.24,0.16],[0.30,0.18],[0.17,0.12],[0.15,0.11],[0.30,0.18],[0.55,0.23],
                 [0.39,0.20],[0.40,0.20],[np.nan,np.nan],[np.nan,np.nan],[np.nan,np.nan]])
best_jetm_err = np.transpose(best_jetm_err)
best_radm_err = np.transpose(best_radm_err)
radm_low_err =   10**(best_radm) - 10**(best_radm-best_radm_err[0])
radm_high_err = 10**(best_radm_err[1]+best_radm) - 10**(best_radm)
jetm_low_err =   10**(best_jetm) - 10**(best_jetm-best_jetm_err[0])
jetm_high_err = 10**(best_jetm_err[1]+best_jetm) - 10**(best_jetm)
plt.rc('font', family='serif')

fig, ax = plt.subplots(1,1, figsize=(7,5))
ax.errorbar(logP_B,10**np.array(best_radm), yerr=[radm_low_err, radm_high_err],
             marker='d',linestyle='none',color ='#018571' ,
             label='Radiative Mode AGN', capsize=3, capthick=1 )
ax.errorbar(logP_B,10**(np.array(best_jetm)), yerr=[jetm_low_err, jetm_high_err],marker='s', linestyle='none',
             color = '#fc8d62', ms = 5, label='Jet Mode AGN', capsize=3, capthick=1 )
logP_5= [24.8, 25.3, 25.8,26.0, 26.5, 27.0 ]
logP_14 = [25.2, 25.7, 26.2, 26.4, 26.9, 27.4]
phi_gps =   10**np.array([-9,-8.15,-8.69,  -9.02, -9.4, -11 ])
phi_gps_err = phi_gps*10**(np.array([-0.4,-0.5,-0.4,-0.3, -0.7,-0.3]))
ax.errorbar(logP_14,phi_gps,yerr=phi_gps_err, marker='*', color ='#535454' , 
              label='GPS', linestyle='none', capsize=3 ,
             capthick=1)

ax.errorbar(L_bins, phi_arr,yerr=phi_err, linestyle='none',marker='o', color='red', 
             label = 'Our Data',capsize=3, capthick=1)


#plt.xscale('log')
ax.set_yscale('log')

labeldict = {'size':16, 'weight':'normal', 'stretch':'normal', 'family':'serif', 
             'style':'normal', 'variant':'normal'}

    


ax.set_xlabel(r'$\rm \log_{10}(L_{1.4}/\,W\,Hz^{-1})$', **labeldict)
ax.set_ylabel(r'$\it\phi\,\,\rm(Mpc^{-3}\,\Delta\log L^{-1})$', **labeldict)

#plt.xlabel('Log'+r'$_{10}(L_{1.4 \rm{GHz}}$' +'/ W Hz'+r'$^{-1})$', **labeldict)
#plt.ylabel('Log('+r'$\Phi$'+' / Mpc'+r'$^{-3} \rm{log L}^{-1})$', **labeldict)

ax.tick_params(axis= 'both', which='major', labelsize=16, length = 7, 
                labeltop = False, top=True, right = True, direction = 'in')
ax.tick_params(axis= 'both', which='minor', labelsize=16, length = 4, labeltop = True, 
                top=True, right = True, direction = 'in')


ax.legend(fontsize=14, frameon=True)
plt.tight_layout()
plt.minorticks_on()
plt.savefig('final_plots/lum_func_v4.pdf')

Vmax_v2 = []
files= glob.glob('Vmax_calc/*.csv')
for fi in files:
    tab = Table.read(fi)
    samps = tab['Sample_Selection']
    