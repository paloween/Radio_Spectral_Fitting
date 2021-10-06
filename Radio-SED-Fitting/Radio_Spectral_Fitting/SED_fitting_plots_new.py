

import sys
sys.path.insert(1, '/Users/pallavipatil/Desktop/VLA/Radio-SED-Fitting/')
import numpy as np
import matplotlib.pyplot as plt
import os

import seaborn as sns
import matplotlib.image as mpimg
import json
import pandas as pd
from matplotlib.patches import Rectangle

from astropy.table import Table,join

import Radio_Models_func
import file_prep_radio_fitting 
import fitting_mcsim_func as mcsim_fit
import SED_cleanup_utils 



loadt = file_prep_radio_fitting.LoadTables()
fprep = file_prep_radio_fitting.FilePrepUtils()
rmfit = Radio_Models_func.RadioModelFit()
sedutils = SED_cleanup_utils.SED_cleanup_utils()

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
#####################
# File with spectral class info:


spdf = pd.read_excel('New_Analysis/Spectral_shape_classf_master.xlsx')
sp_class = Table.from_pandas(spdf)
sp_class.add_index('Source_name')

final_unq_jvla = Table.read('New_Analysis/final_unq_JVLASP_vf_morph_regmod.csv')
final_unq_jvla.add_index('Source_name')

fit_results_tab = Table.read('New_Analysis/Fit_results_new_tab.csv')
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

al_calcs = Table.read('New_Analysis/Radio_color_color_SP_class_info.csv')
al_calcs.add_index('Source_name')

tab = pd.DataFrame({'survey' : ['VLA-A/B', '', 'VLASS',  'NVSS', 'FIRST', 'SUMSS', 'WENSS', 
          'GLEAM', 'TGSS', 'VLSSr', 'VLITE'], 
'cen_freq': [10, 10, 3, 1.4, 1.4, 0.843, 0.325, 0.200, 0.147, 0.074, 0.364],
'sensitivity' : [0.013, 0.013, 0.15, 0.45, 0.15, 1.1, 3., 8, 3.5, 100, 5],
'res': [0.2, 0.6, 2.5, 45, 5, 45, 54, 100, 25, 75, 10], 
'bw': [2, 2, 1, 0.03, 0.03,0.002, 0.003, [[0.126], [0.031]],0.0085,  0.002, 0.032 ]})
df = Table.from_pandas(tab)
df.add_index('survey')

vlite_sources = {'J1238+52':17,  'J0823-06':16.2, 'J0642-27':20 }

fit_results_ext_tab = Table.read('Fit_results_new_VCSS_RACS_inband_fscale.csv')
fit_results_ext_tab.add_index('source')

vcss_cat = Table.read('../Radio-SED-Fitting/Datasets/Pallavi_vcss_mosaics_source_matches_5arcsec.fits')
vcss_cat.add_index('Source_name')

vcss_cat_sources = []
for row in vcss_cat:
    vcss_cat_sources.append(row['Source_name'].rstrip())
     
racs_tab = Table.read('../Radio-SED-Fitting/Datasets/RACS-DR1_Xmatch_WISE.fits')
if 'Source_name_1' in racs_tab.colnames:
    racs_tab.rename_column('Source_name_1', 'Source_name')
racs_tab.add_index('Source_name')


atgcat.add_index('WISEname')
atscat.add_index('WISEname')
sp_class.add_index('Source_name')

myfdf = pd.read_csv('Tab4_Flux_measurements_spectra_fscale.csv')

plot_flag = 'all'
plt.ioff()

models = ['PL', 'CPL', 'EFFA', 'IFFA', 'SSA', 'GCA']
colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
#rgb = sns.xkcd_palette(colors)
rgb = sns.color_palette("Set1", 10)
sns.set(style="ticks", palette="pastel", color_codes=True)
plt.rcParams.update({'lines.markeredgewidth': 1})
lincols = ['#66c2a5', '#fc8d62','#8da0cb','#e78ac3','#a6d854']

plot_conf = False
plot_parreg = False
plotrcolor = True
plotmc = False

plotgcv= False
diagmode = True


for ii, row in enumerate(final_unq_jvla):
    source = row['Source_name']
    if source != 'J0457-23':       
        wise_row = wise_phot.loc['Source_name',source]
        wisen = wise_row['WISEname']
        sp_row = sp_class.loc[source]
        morph = sp_row['Xband_morph']
        lomorph = sp_row['Visual_morph']
        spshape = sp_row['Modified_Sp_Shape']
        mc_par = par_fits[source]
        scat = wise_row
        glmcat = atgcat.loc[wisen]
    
        jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
        jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
         
        freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
        OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,jvla_BX)
        freq_arr, flux_arr, eflux_arr, labels = sedutils.modify_data_plot(source,freq_arr,
                            flux_arr, eflux_arr, labels,sp_row, wise_row ) 
        
        #get VCSS data
        #freq_arr, flux_arr, eflux_arr, labels = get_VCSS_data(source, wise_row,
        #                     sp_row, labels, freq_arr, flux_arr, eflux_arr)
        freq_arr, flux_arr, eflux_arr, labels = sedutils.VCSS_data_cleanup(source, wise_row,
                            sp_row, labels, freq_arr, flux_arr, eflux_arr)
    
        flux_arr, eflux_arr, labels = sedutils.mod_data_flux_density_scale(flux_arr, eflux_arr, 
                                                              labels)
 
        
        if 'VLASS' in labels and sp_row['VLASS_Limit']!='U':
            ind = labels.index('VLASS')
            vflux = flux_arr[ind]
            eflux = eflux_arr[ind]
            if eflux<0.2*vflux:
                eflux_arr[ind] = np.sqrt(eflux**2 + (vflux/5)**2)        
        
        name = scat['WISEname']  #Name of the object.
        msk = [freq_arr>0]
        nu_arr = freq_arr[tuple(msk)]
        s_nu= flux_arr[tuple(msk)]
        e_s_nu = eflux_arr[tuple(msk)]
        survey_faint = ['TGSS', 'WENSS', 'GLEAM', 'SUMSS', 'VLSSr']   
        freq_faint = [0.15, 0.325, 0.200, 0.843, 0.074]
        for survey, freq in zip(survey_faint, freq_faint):
            flag = sp_row[survey+'_Limit']
            if flag == 'T' or flag=='F':
                if survey not in labels:
                    f_flux, f_eflux = mcsim_fit.get_faint_detection(source, survey,10)
                    if f_flux>0:
                        labels.append(survey)
                        nu_arr = np.append(nu_arr, freq)
                        s_nu = np.append(s_nu, f_flux)
                        e_s_nu = np.append(e_s_nu, f_eflux)

        nreg = max(1, len(jvla_AX), len(jvla_BX))
        if source == 'J2325-04': nreg = 1
        # Adding individual points for multi-coomponent sources.
        config = row['config']
        mystr = '{0}\nz: {1}\nX Morph: {2}\nLo Morph: {3}\nSED: {4}\nconfig: {5}\n'.format(
                               source,sp_row['redshift'], morph, lomorph, spshape, config )
        mystr = spshape
        if sp_row['redshift'] >0:
            mystr = mystr+', z: {0}'.format(sp_row['redshift'])

        fitRes = fit_results_ext_tab.loc[source]
        par_pl = fitRes['PL_s0', 'PL_al']
        par_cpl = fitRes['CPL_s0', 'CPL_al', 'CPL_nup']
        par_effa = fitRes['EFFA_s0', 'EFFA_al', 'EFFA_nup']
        par_iffa = fitRes['IFFA_s0', 'IFFA_al', 'IFFA_nup']
        par_ssa = fitRes['SSA_s0', 'SSA_al', 'SSA_nup']
        perr_pl = fitRes['perr_pl_s0','perr_pl_al']
        perr_cpl = fitRes['perr_cpl_s0','perr_cpl_al','perr_cpl_nup']
        
        
        mcpar_cpl = [mc_par['s0_cpl']['Median'], mc_par['al_cpl']['Median'], mc_par['q_cpl']['Median']]
        mcpar_cpl_err = [mc_par['s0_cpl']['Std dev'], mc_par['al_cpl']['Std dev'], mc_par['q_cpl']['Std dev']]
        
        nu_mod = np.logspace(-3, 3.0, 5000)
                
        s1 = rmfit.power_law(nu_mod, *par_pl)
        s2 = rmfit.CPL(nu_mod, *par_cpl)
        s3 = rmfit.EFFA_func(nu_mod, *par_effa )
        s4 = rmfit.IFFA_func(nu_mod, *par_iffa)
        s5 = rmfit.SSA_func(nu_mod, *par_ssa)
        s6 = rmfit.CPL(nu_mod, *mcpar_cpl)                      

        if plotgcv:
            par_gcv = fitRes['GCV_s0', 'GCV_al1', 'GCV_al2', 'GCV_nup']
            s7 = rmfit.gen_curve_func(nu_mod, *par_gcv)  
            
        cpl_al = par_cpl['CPL_al']
        cpl_q = par_cpl['CPL_nup']
        cpl_nup = np.exp(-cpl_al/(2*cpl_q))
        cpl_nup_err = abs((-cpl_nup/(2*cpl_q))*(perr_cpl[1] + (cpl_al*perr_cpl[2]/cpl_q)))
 
        par_str = r' $\alpha_{PL}:$'+' {0:.2f},'.format(par_pl['PL_al']) 
        par_str = par_str + r' $\alpha_{CV}:$'+' {0:.2f}'.format(par_cpl['CPL_al']) +r'$\pm$'+' {0:.2f}\n'.format(perr_cpl[1])
        par_str = par_str+r' $\nu^{p}_{CV}:$'+' {0:.1e}'.format(cpl_nup)+r'$\pm$'+' {0:.1e} GHz\n'.format(cpl_nup_err)
        par_str = par_str+r' $q_{CV}:$'+' {0:.1f}'.format(cpl_q)+r'$\pm$'+' {0:.1f}'.format(perr_cpl[2])
        par_str = par_str+r' $S_{0,CV}:$'+' {0:.1f}'.format(par_cpl[0])+r'$\pm$'+' {0:.1f}\n'.format(perr_cpl[0])
        par_str = par_str+r'$\chi^{2}_{r, PL}:$'+'{0:.1f}, '.format(fitRes['rchisq_pl'])
        par_str = par_str+r'$\chi^{2}_{r, CPL}:$'+'{0:.1f}\n'.format(fitRes['rchisq_cpl'])

        
        if plotmc:
            pars, cpl_fit_val, tcarlo  = mcsim_fit.mcsim_fitting(wise_row, 'CPL', atgcat, vla_ax_grp, 
                                                            vla_bx_grp, 1000)
            cpl_fit =  dict.fromkeys(['s0_cpl', 'al_cpl', 'q_cpl', 'speak', 'nupeak' ])
            for key, prow in zip(cpl_fit.keys(), cpl_fit_val):
                pardict = {} 
                for k, v in zip(pars, prow):
                    pardict[k] = v
                cpl_fit[key] = pardict            
            med_pars = [np.median(tcarlo[0][:]),np.median(tcarlo[1][:]),np.median(tcarlo[2][:])]
            med_epar = [np.std(tcarlo[0][:]),np.std(tcarlo[1][:]),np.std(tcarlo[2][:])]

            par_str = par_str + r' $\alpha_{CV,M}:$'+' {0:.2f}'.format(med_pars[1])+r'$\pm$'+' {0:.2f},'.format(med_epar[1])
            par_str = par_str+r' $q_{CV,M}:$'+' {0:.1f}\n'.format(med_pars[2])+r'$\pm$'+' {0:.2f},\n'.format(med_epar[2])
            par_str = par_str+r' $S_{0,CV,M}:$'+' {0:.1f}'.format(med_pars[0])+r'$\pm$'+' {0:.2f}'.format(med_epar[0])

        fig, ax= plt.subplots(1,1,figsize= (6.52, 8.34))

        if 'AX' in labels and len(alpha_AX)>0:
            if np.any(jvla_AX['snr']>50) and alpha_AX[1]>-999:
                mind = labels.index('AX')
                FAX_10GHZ = flux_arr[mind]
                EAX_10GHZ = eflux_arr[mind]
                sedutils.plot_bowtie(10,alpha_AX, FAX_10GHZ, EAX_10GHZ,ax, zorder=9)
        if 'BX' in labels and len(alpha_BX)>0:
            if np.any(jvla_BX['snr']>50) and alpha_BX[1]>-999:
                mind = labels.index('BX')
                FBX_10GHZ = flux_arr[mind]
                EBX_10GHZ = eflux_arr[mind]
                sedutils.plot_bowtie(10,alpha_BX, FBX_10GHZ, EBX_10GHZ,ax, zorder=9)
            
                
        if sp_row['TGSS_Limit'] =='U':
            noise = wise_row['TGSS_Limit']
            sedutils.plot_upper_limit(ax, 0.15, 4*noise, 3*noise, 9)
            
        if sp_row['VLSSr_Limit'] =='U':
            labels.append('VLSSr')
            noise = wise_row['VLSSr_Limit']
            sedutils.plot_upper_limit(ax, 0.074, 4*noise, 3*noise, 7)

        if sp_row['WENSS_Limit'] =='U':            
            noise = wise_row['WENSS_Limit']
            if noise>0:
                labels.append('WENSS')
                sedutils.plot_upper_limit(ax, 0.325, 4*noise, 3*noise, 7)
                
        if sp_row['GLEAM_Limit'] =='U':            
            noise = wise_row['GLEAM_Limit']
            if noise>0:
                labels.append('GLEAM')
                sedutils.plot_upper_limit(ax, 0.200, 3*noise, 2*noise, 7)
        #adding VCSS upper limit:
        
        
                
       # vlite_det = False
       # if 'VLITE' in labels:
       #     my_ind = labels.index('VLITE')
       #     vlite_det  = True
       #     vlite_freq = nu_arr[my_ind]
       #     vlite_flux = s_nu[my_ind]
       #     vlite_eflux= e_s_nu[my_ind]
       #     nu_arr = np.delete(nu_arr, my_ind)
       #     s_nu = np.delete(s_nu, my_ind)
       #     e_s_nu = np.delete(e_s_nu, my_ind)
                    
        ax.errorbar(nu_arr, s_nu, yerr=e_s_nu, marker ='o', 
                    linestyle='none', capsize=3, markersize = 4, color='k', 
                    elinewidth = 1, zorder = 10)
        
        #if source in list(racs_tab['Source_name']):
        #    racs_row = racs_tab.loc[source]
        #    ax.errorbar(0.887, racs_row['flux_int'], marker ='*', 
        #            linestyle='none',  markersize = 8, color='tab:pink', 
        #            elinewidth = 1, zorder = 10)
        #    ax.axvspan(0.743, 1.031,facecolor='tab:pink', alpha=0.1, zorder=1 )
        #    ax.annotate('RACS', xy=(.46,.2), xytext=(.46,0.2), color='tab:pink', 
        #            fontsize=12, fontweight='medium',xycoords='axes fraction' )
            
            
        ax.text(0.03, 0.9,mystr, fontsize=14, bbox=dict(facecolor='white', alpha=0.5, 
                                          boxstyle='round', ec = 'gray'), 
                                transform=ax.transAxes, zorder = 11)
            
           
        # Adjusting the labels, axes scales and range.  
        ax.set_title(source,  fontsize=16, fontweight='medium')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(0.01*min(flux_arr), 10*max(flux_arr))
        ax.set_xlim(0.03, 20)
        ax.set_xlabel(r'log ($\nu / GHz)$', fontdict=dict(fontsize=14, ))
        ax.set_ylabel(r'log ($S_{\nu}/mJy)$', fontdict=dict(fontsize=14))
        ax.tick_params(axis = 'both', which = 'minor', direction='in', 
                       length=4, top=True, right=True)
        ax.tick_params(axis = 'both', which = 'major', direction='in',
                       length=9, top=True, right=True, labelsize=14)
    
        
            
    
        
        if plotmc:
            for p1, p2, p3 in zip(tcarlo[0][:],tcarlo[1][:],tcarlo[2][:] ):
                s7 = rmfit.CPL(nu_mod, *[p1, p2, p3])   
                ax.plot(nu_mod, s7, linestyle = '-', color='gray', 
                            lw=0.5, alpha = 0.1, zorder = 2)
            s6 = rmfit.CPL(nu_mod, *med_pars) 
            l6, = ax.plot(nu_mod, s6, linestyle = '-', label='CV-MC', color=lincols[2], 
                      lw=2, zorder = 3)

        if plotgcv:
            l7, = ax.plot(nu_mod, s7, linestyle = '-', label='GCV', color=lincols[2], 
                      lw=2, zorder = 3)
        
        fitRes = fit_results_ext_tab.loc[source]
        par_pl = fitRes['PL_s0', 'PL_al']
        par_cpl = fitRes['CPL_s0', 'CPL_al', 'CPL_nup']
        
        nu_mod = np.logspace(-3, 3.0, 5000)
        
        s1 = rmfit.power_law(nu_mod, *par_pl)
        s2 = rmfit.CPL(nu_mod, *par_cpl)
        lincols = ['#66c2a5', '#fc8d62','#8da0cb','#e78ac3','#a6d854']
                   
        if diagmode:
            l1, = ax.plot(nu_mod, s1, linestyle = 'dotted', label='PL', color=lincols[0], 
                      lw=1, zorder = 1)     
            l2, = ax.plot(nu_mod, s2, linestyle = 'dotted', label='CV', color=lincols[1], 
                      lw=1, zorder = 1)
            ax.text(0.48, 0.77,par_str, fontsize=9, bbox=dict(facecolor='white', alpha=0.5, 
                                      boxstyle='round', ec = 'gray'), 
                            transform=ax.transAxes, zorder=12)
     
        else:
            if spshape in ['PL', 'PL:', 'I', 'F']:
                l1, = ax.plot(nu_mod, s1, linestyle = 'dotted', label='PL', color=lincols[0], 
                          lw=1, zorder = 1)       
            else:
                l2, = ax.plot(nu_mod, s2, linestyle = 'dotted', label='CV', color=lincols[1], 
                          lw=1, zorder = 1)
        
           
        ax.axvspan(1.37, 1.43,facecolor=rgb[0], alpha=0.2, zorder=1 )
        ax.axvspan(8, 12,facecolor=rgb[3], alpha=0.2, zorder=1 )
        ax.axvspan(.139, .156,facecolor=rgb[4], alpha=0.2, zorder=1 )
        ax.axvspan(2, 4,facecolor=rgb[8], alpha=0.2, zorder=1 )
        
    
        ax.annotate('NVSS', xy=(.51,.2), xytext=(.55,0.2), color=rgb[0], 
                    fontsize=12, fontweight='medium',xycoords='axes fraction' )
        ax.annotate('TGSS', xy=(.13,.2), xytext=(.20,0.2), color=rgb[4], 
                    fontsize=12, fontweight='medium',xycoords='axes fraction' )    
        ax.annotate('VLA-X', xy=(.83,.2), xytext=(.85,.2), color=rgb[3], 
                    fontsize=12, fontweight='medium',xycoords='axes fraction' )
        if (wisen =='040403.61-243600.1') or (wisen == '061200.23-062209.1'):
            ax.annotate('VLASS', xy=(.64,.1), xytext=(.67,.1), color=rgb[8], 
                    fontsize=12, fontweight='medium',xycoords='axes fraction' )
        else:
            ax.annotate('VLASS', xy=(.64,.12), xytext=(.65,.12), color=rgb[8], 
                    fontsize=12, fontweight='medium',xycoords='axes fraction' )
        
        if 'VLITE' in labels:
            mycolor = 'tab:red'
            #ax.errorbar(vlite_freq, vlite_flux, yerr=vlite_eflux, marker ='*', 
            #    linestyle='none', capsize=3, markersize = 4, color=mycolor, 
            #    elinewidth = 1, zorder = 10)
    
            ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
            ax.annotate('VCSS', xy=(.31,.2), xytext=(.33,0.2), color=mycolor, 
                fontsize=12, fontweight='light',xycoords='axes fraction' )   
        
        if ii not in [103,139]:
            if 'WENSS' in labels or 'TEXAS' in labels or 'WISH' in labels:
                ax.axvspan(.322, .328,facecolor=rgb[2], alpha=0.2 )
            #ax.annotate('VLITE', xy=(.28,.2), xytext=(.28,0.2), color=rgb[2], 
            #        fontsize=12, fontweight='medium',xycoords='axes fraction' )
        else:
            ax.axvspan(.35, .38,facecolor=rgb[2], alpha=0.2, zorder=1 )
       
        if 'GLEAM' in labels: 
            ax.axvspan(0.17, 0.23,facecolor=rgb[6], alpha=0.1, zorder=1 )
            ax.annotate('GLEAM', xy=(.16,.12), xytext=(.23,0.12), color=rgb[6], 
                    fontsize=12, fontweight='light',xycoords='axes fraction' ) 
        if 'WENSS' in labels:
            ax.annotate('WENSS', xy=(.31,.12), xytext=(.31,0.12), color=rgb[2], 
                    fontsize=12, fontweight='light',xycoords='axes fraction' ) 
        if 'VLSSr' in labels:
            ax.axvspan(.073, .075,facecolor=rgb[1], alpha=0.1, zorder=1)
            ax.annotate('VLSSr', xy=(.05,.12), xytext=(.06,0.12), color=rgb[1], 
                    fontsize=12, fontweight='light',xycoords='axes fraction' )
        if source =='J1651+34':
            ax.axvspan(3.85, 5.85,facecolor='gray', alpha=0.1, zorder=1 )
            ax.annotate('GB6', xy=(.75,.2), xytext=(.75,0.2), color='gray', 
                    fontsize=12, fontweight='light',xycoords='axes fraction' )
        if 'RACS' in labels:
            ax.axvspan(0.743, 1.031,facecolor=rgb[8], alpha=0.1, zorder=1 )
            ax.annotate('RACS', xy=(.46,.12), xytext=(.46,0.12), color=rgb[8], 
                        fontsize=12, fontweight='medium',xycoords='axes fraction' )
           
        if 'Texas' in labels:
            ax.annotate('TEXAS', xy=(.32,.2), xytext=(.33,0.2), color=rgb[2], 
                    fontsize=12, fontweight='medium',xycoords='axes fraction' )
            '''
            else:    
                if sp_row['VLITE_Limit']=='F':               
                    mycolor= 'tab:cyan'
                    ax.errorbar(vlite_freq, vlite_flux, yerr=vlite_eflux, marker ='*', 
                        linestyle='none', capsize=3, markersize = 4, color=mycolor, 
                        elinewidth = 1, zorder = 10)
                    ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
                    ax.annotate('VCSS', xy=(.31,.12), xytext=(.33,0.12), color=mycolor, 
                        fontsize=12, fontweight='light',xycoords='axes fraction' )
    
                 
                elif wise_row['VLITE_Limit']>0:
                    mycolor = 'tab:cyan'
                    plot_upper_limit(0.338, 2.5*wise_row['VLITE_Limit'], 
                            1.5*wise_row['VLITE_Limit'], 7,color=mycolor)
                    ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
                    ax.annotate('VCSS', xy=(.31,.12), xytext=(.33,0.12), color=mycolor, 
                        fontsize=12, fontweight='light',xycoords='axes fraction' )
        '''
        '''
        if plot_conf:
           try:
               pl_upp, pl_low, pl_err = calc_ci(par_pl, perr_pl, nu_mod, rmfit.power_law)
               ax.fill_between(nu_mod, pl_low, pl_upp, fc =lincols[0], alpha = 0.3 )
           except:
                continue
           try:
               cpl_upp, cpl_low, cpl_err = calc_ci(par_cpl, perr_cpl, nu_mod,
                                            rmfit.CPL)
               ax.fill_between(nu_mod, cpl_low, cpl_upp, fc =lincols[1], alpha = 0.3 )
           except: continue
        '''
        
        if plot_parreg and abs(med_pars[2])>0.1:
            lpert = 15.9
            hpert = 84.1

            x = cpl_fit['nupeak']['Median']-cpl_fit['nupeak']['e1Low']
            y = cpl_fit['speak']['Median']-cpl_fit['speak']['e1Low']
            dx = cpl_fit['nupeak']['e1Low']+cpl_fit['nupeak']['e1High']
            dy = cpl_fit['speak']['e1Low']+cpl_fit['speak']['e1High']
            
            #x = cpl_fit['speak']['Median'] - cpl_fit['speak']['e1Low']
            #x = cpl_fit['speak']['Median'] - cpl_fit['speak']['e1Low']

            rect_cpl = Rectangle((x,y),dx, dy, fc = lincols[1], 
                               label='CV-MC', alpha = 0.5,ec='white', linewidth=2.0 )            
            try:
                ax.add_patch(rect_cpl)
            except: 
                pass
        

        
        #ax.annotate(name,xy=(.3,1.1), xytext=(.3,1.1), color='k', 
        #         fontsize=12, fontweight='light' ,xycoords='axes fraction' )
        ax.legend(loc='upper right')
        plt.tight_layout()
        plt.subplots_adjust(bottom = 0.45, right = 0.99)
        ax3 = plt.axes([0.02, 0.01, 0.97, 0.35])
        ax3.set_xticks([])
        ax3.set_yticks([])
        imgfile = 'figures/'+ source.strip()+'_comb_v6.png'
        if os.path.exists(imgfile):
            img = mpimg.imread(imgfile)
            ax3.imshow(img, aspect='auto')
            ax3.set_xticks([])
            ax3.set_yticks([])
            
            
        # Showing Radio Color-Color Plot:
        if plotrcolor:
            axins = ax.inset_axes([0.07, 0.28, 0.2, 0.2])
            alrow = al_calcs.loc[source]
            x = alrow['al_l_TGSS_NVSS']
            y = alrow['al_h_NVSS_VLAX']
            xerr = alrow['e_al_l_TGSS_NVSS']
            yerr = alrow['e_al_h_NVSS_VLAX']
            lim = alrow['l_al_l_TGSS_NVSS']
            if lim.lower() == 'true': 
                lim = True
            else:  lim = False
            axins.errorbar(x, y, xerr=xerr, yerr=yerr, marker='o', 
                        linestyle='none', capsize=2, ms=4, c='tab:blue', 
                        xlolims = lim)
             
            axins.axhline(0, ls='-', c='gray', lw=1.0, zorder=1)
            axins.axvline(0, ls='-', c='gray', lw=1.0, zorder=1)
            axins.axhline(-0.5, ls='--', c='purple', lw=1.0, alpha=0.8, zorder=1)
            axins.axvline(0.1, ls='--', c='purple', lw=1.0, alpha=0.8, zorder=1)
            
            
            x = np.linspace(-3.0, 3.0, 50)
            axins.plot(x, x, ls='-', c='gray', lw=1.0, zorder=1)
            axins.set_xlim(-1.5, 1.8)
            axins.set_ylim(-2.3, 0.6)
            axins.set_xlabel(r'$\alpha^{1.4}_{0.15}$', fontsize=10, 
                             labelpad = 0.01)
            axins.set_ylabel(r'$\alpha^{10}_{1.4}$', fontsize=10, 
                             labelpad = 0.01)
            axins.set_xlim(-2.0, 2.0)
            axins.set_ylim(-3.0, 2.0)
            axins.set_xticks([-2,0,2])
            axins.set_yticks([-2,0,2])
            axins.tick_params(which='major', labelsize=10, direction='in',length=8, 
                        top=True, right=True, labelbottom=False, labelleft=False,
                        labeltop=True,labelright=True)
            
        plt.savefig('RadioSED_plots/New_Fit_Color_Plots/'+name+'_rsed_20.pdf',
                    dpi = 100)
        
        plt.close()
        
        
        
        
        
'''
SED cutouts
'''        
  
inband_fit = Table.read('Fit_results_new_VCSS_RACS_inband_fscale.csv')      
inband_fit.add_index('source')
nu_mod = np.logspace(-1.3,1.3, 100 )
nrows = 5
ncols = 3
filecounter = 0
myfilecounter = []
subplotcounter = 0
mysubplotcounter = []
plotarr = [np.nan]*(nrows*ncols)
npage = 0
nplots = nrows*ncols
plt.ioff()
rgb = sns.color_palette("Set1", 10)

sources_ib =['J0417-28', 'J0525-36', 'J2251+01']

for ii, row in enumerate(final_unq_jvla): 
    source = row['Source_name']
    wise_row = wise_phot.loc['Source_name',source]
    wisen = wise_row['WISEname']
    al_row = al_calcs.loc[source]
    sp_row = sp_class.loc[source]
    spshape = sp_row['Modified_Sp_Shape']
    morph = al_row['Xband_morph']
    
    curved_cats = ['CV', 'CV:', 'PK', 'PK:']
    model_name = 'PL'
    if np.isin(spshape, curved_cats):
        model_name = 'CPL'
    
    
    scat = wise_row
    glmcat = atgcat.loc[wisen]
    
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
    
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
    OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,jvla_BX)
    freq_arr, flux_arr, eflux_arr, labels = sedutils.modify_data_plot(source,freq_arr,
                        flux_arr, eflux_arr, labels,sp_row, wise_row )        
    
    freq_arr, flux_arr, eflux_arr, labels = sedutils.VCSS_data_cleanup(source, wise_row,
                            sp_row, labels, freq_arr, flux_arr, eflux_arr)
    
    if 'VLASS' in labels and sp_row['VLASS_Limit']!='U':
        ind = labels.index('VLASS')
        vflux = flux_arr[ind]
        eflux = eflux_arr[ind]
        if eflux<0.2*vflux:
            eflux_arr[ind] = np.sqrt(eflux**2 + (vflux/5)**2)        
    
    
    name = scat['WISEname']  #Name of the object.
    msk = [freq_arr>0]
    nu_arr = freq_arr[tuple(msk)]
    s_nu= flux_arr[tuple(msk)]
    e_s_nu = eflux_arr[tuple(msk)]
    survey_faint = ['TGSS', 'WENSS', 'GLEAM', 'SUMSS', 'VLSSr']   
    freq_faint = [0.15, 0.325, 0.200, 0.843, 0.074]
    for survey, freq in zip(survey_faint, freq_faint):
        flag = sp_row[survey+'_Limit']
        if flag == 'T' or flag=='F':
            if survey not in labels:
                f_flux, f_eflux = mcsim_fit.get_faint_detection(source, survey,10)
                if f_flux>0:
                    labels.append(survey)
                    nu_arr = np.append(nu_arr, freq)
                    s_nu = np.append(s_nu, f_flux)
                    e_s_nu = np.append(e_s_nu, f_eflux)
    
    nreg = max(1, len(jvla_AX), len(jvla_BX))
    if source == 'J2325-04': nreg = 1
    
    ymin = min(flux_arr)
    ymax = max(flux_arr)
    
    axind = ii%nplots
    axr = int(axind/ncols)
    axc = axind - axr*ncols
    print(axind, axr, axc)
    if axind ==0:
        fig,axs = plt.subplots(nrows,ncols, num=npage, figsize = (15,19))     
    ax = axs[axr, axc]
    
    
    if 'AX' in labels and len(alpha_AX)>0:
        if np.any(jvla_AX['snr']>50) and alpha_AX[1]>-999:
            mind = labels.index('AX')
            FAX_10GHZ = flux_arr[mind]
            EAX_10GHZ = eflux_arr[mind]
            sedutils.plot_bowtie(10,alpha_AX, FAX_10GHZ, EAX_10GHZ,ax, zorder=9)
    if 'BX' in labels and len(alpha_BX)>0:
        if np.any(jvla_BX['snr']>50) and alpha_BX[1]>-999:
            mind = labels.index('BX')
            FBX_10GHZ = flux_arr[mind]
            EBX_10GHZ = eflux_arr[mind]
            sedutils.plot_bowtie(10,alpha_BX, FBX_10GHZ, EBX_10GHZ,ax, zorder=9)           
    
    if sp_row['TGSS_Limit'] =='U':
        noise = wise_row['TGSS_Limit']
        sedutils.plot_upper_limit(ax, 0.15, 2*noise, noise, 9)
        ymin = min(ymin, noise)
        ymax = max(ymax, 3*noise)
    
    if sp_row['VLSSr_Limit'] =='U':
        labels.append('VLSSr')
        noise = wise_row['VLSSr_Limit']
        sedutils.plot_upper_limit(ax, 0.074, 2*noise, noise, 7)
        ymin = min(ymin, noise)
        ymax = max(ymax, 3*noise)
    
    if sp_row['WENSS_Limit'] =='U':            
        noise = wise_row['WENSS_Limit']
        if noise>0:
            labels.append('WENSS')
            sedutils.plot_upper_limit(ax, 0.325, 2*noise, noise, 7)
            ymin = min(ymin, noise)
            ymax = max(ymax, 3*noise)
    
    if sp_row['GLEAM_Limit'] =='U':            
        noise = wise_row['GLEAM_Limit']
        if noise>0:
            labels.append('GLEAM')
            sedutils.plot_upper_limit(ax, 0.200, 2*noise, noise, 7)
            ymin = min(ymin, noise)
            ymax = max(ymax, 3*noise)
    
    ax.errorbar(nu_arr, s_nu, yerr=e_s_nu, marker ='o', 
                linestyle='none', capsize=3, markersize = 4, color='k', 
                elinewidth = 1, zorder = 10)            
    
    
    ax.set_title(source,  fontsize=16, fontweight='medium')
    ax.set_xscale('log')
    ax.set_yscale('log')
    if spshape=='PL':
        ax.set_ylim(ymin/2, 3*ymax)
    else:
        ax.set_ylim(ymin/2, 3*ymax)
    ax.set_xlim(0.03, 20)
    ax.set_xlabel(r'log ($\nu / GHz)$', fontdict=dict(fontsize=14, ))
    ax.set_ylabel(r'log ($S_{\nu}/mJy)$', fontdict=dict(fontsize=14))
    ax.tick_params(axis = 'both', which = 'minor', direction='in', 
                   length=4, top=True, right=True)
    ax.tick_params(axis = 'both', which = 'major', direction='in',
                   length=9, top=True, right=True, labelsize=14)
    
    
    ax.text(0.9,0.9, spshape, 
            transform=ax.transAxes, color='k', 
            fontsize = 14,bbox = dict(boxstyle='round',fc='w', ec='k'))
    ax.text(0.8,0.9, morph, 
            transform=ax.transAxes, color='k', 
            fontsize = 14,bbox = dict(boxstyle='round',fc='w', ec='k'))
    
    
    fitRes = fit_results_ext_tab.loc[source]
    if source in sources_ib:
        fitRes = inband_fit.loc[source]
    par_pl = fitRes['PL_s0', 'PL_al']
    par_cpl = fitRes['CPL_s0', 'CPL_al', 'CPL_nup']
    
    nu_mod = np.logspace(-3, 3.0, 5000)
    
    s1 = rmfit.power_law(nu_mod, *par_pl)
    s2 = rmfit.CPL(nu_mod, *par_cpl)
    lincols = ['#66c2a5', '#fc8d62','#8da0cb','#e78ac3','#a6d854']
    
    if spshape in ['PL', 'PL:', 'I', 'F']:
        l1, = ax.plot(nu_mod, s1, linestyle = 'dotted', label='PL', color=lincols[0], 
                  lw=1, zorder = 1)       
    else:
        l2, = ax.plot(nu_mod, s2, linestyle = 'dotted', label='CV', color=lincols[1], 
                  lw=1, zorder = 1)
    
    z = al_row['redshift']
    if z>0:
        ax.text(0.8, 0.8,r'$z$'+'={0}'.format(z), fontsize=14, bbox=dict(facecolor='white', alpha=0.5, 
                                  boxstyle='round', ec = 'gray'), 
                        transform=ax.transAxes, zorder = 11)
    
    ax.axvspan(1.37, 1.43,facecolor=rgb[0], alpha=0.2, zorder=1 )
    ax.axvspan(8, 12,facecolor=rgb[3], alpha=0.2, zorder=1 )
    ax.axvspan(.139, .156,facecolor=rgb[4], alpha=0.2, zorder=1 )
    ax.axvspan(2, 4,facecolor=rgb[8], alpha=0.2, zorder=1 )
    

    ax.annotate('NVSS', xy=(.51,.2), xytext=(.55,0.2), color=rgb[0], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    ax.annotate('TGSS', xy=(.13,.2), xytext=(.20,0.2), color=rgb[4], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )    
    ax.annotate('VLA-X', xy=(.83,.2), xytext=(.85,.2), color=rgb[3], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    if (wisen =='040403.61-243600.1') or (wisen == '061200.23-062209.1'):
        ax.annotate('VLASS', xy=(.64,.1), xytext=(.67,.1), color=rgb[8], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    else:
        ax.annotate('VLASS', xy=(.64,.12), xytext=(.65,.12), color=rgb[8], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    
    if 'VLITE' in labels:
        mycolor = 'tab:red'
        #ax.errorbar(vlite_freq, vlite_flux, yerr=vlite_eflux, marker ='*', 
        #    linestyle='none', capsize=3, markersize = 4, color=mycolor, 
        #    elinewidth = 1, zorder = 10)

        ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
        ax.annotate('VCSS', xy=(.31,.2), xytext=(.33,0.2), color=mycolor, 
            fontsize=12, fontweight='light',xycoords='axes fraction' )   
    
    if ii not in [103,139]:
        if 'WENSS' in labels or 'TEXAS' in labels or 'WISH' in labels:
            ax.axvspan(.322, .328,facecolor=rgb[2], alpha=0.2 )
        #ax.annotate('VLITE', xy=(.28,.2), xytext=(.28,0.2), color=rgb[2], 
        #        fontsize=12, fontweight='medium',xycoords='axes fraction' )
    else:
        ax.axvspan(.35, .38,facecolor=rgb[2], alpha=0.2, zorder=1 )
   
    if 'GLEAM' in labels: 
        ax.axvspan(0.17, 0.23,facecolor=rgb[6], alpha=0.1, zorder=1 )
        ax.annotate('GLEAM', xy=(.16,.12), xytext=(.23,0.12), color=rgb[6], 
                fontsize=12, fontweight='light',xycoords='axes fraction' ) 
    if 'WENSS' in labels:
        ax.annotate('WENSS', xy=(.31,.12), xytext=(.31,0.12), color=rgb[2], 
                fontsize=12, fontweight='light',xycoords='axes fraction' ) 
    if 'VLSSr' in labels:
        ax.axvspan(.073, .075,facecolor=rgb[1], alpha=0.1, zorder=1)
        ax.annotate('VLSSr', xy=(.05,.12), xytext=(.06,0.12), color=rgb[1], 
                fontsize=12, fontweight='light',xycoords='axes fraction' )
    if source =='J1651+34':
        ax.axvspan(3.85, 5.85,facecolor='gray', alpha=0.1, zorder=1 )
        ax.annotate('GB6', xy=(.75,.2), xytext=(.75,0.2), color='gray', 
                fontsize=12, fontweight='light',xycoords='axes fraction' )
    if 'RACS' in labels:
        ax.axvspan(0.743, 1.031,facecolor=rgb[8], alpha=0.1, zorder=1 )
        ax.annotate('RACS', xy=(.46,.12), xytext=(.46,0.12), color=rgb[8], 
                    fontsize=12, fontweight='medium',xycoords='axes fraction' )
       
    if 'Texas' in labels:
        ax.annotate('TEXAS', xy=(.32,.2), xytext=(.33,0.2), color=rgb[2], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    '''
        else:    
            if sp_row['VLITE_Limit']=='F':               
                mycolor= 'tab:cyan'
                ax.errorbar(vlite_freq, vlite_flux, yerr=vlite_eflux, marker ='*', 
                    linestyle='none', capsize=3, markersize = 4, color=mycolor, 
                    elinewidth = 1, zorder = 10)
                ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
                ax.annotate('VCSS', xy=(.31,.12), xytext=(.33,0.12), color=mycolor, 
                    fontsize=12, fontweight='light',xycoords='axes fraction' )

             
            elif wise_row['VLITE_Limit']>0:
                mycolor = 'tab:cyan'
                plot_upper_limit(0.338, 2.5*wise_row['VLITE_Limit'], 
                        1.5*wise_row['VLITE_Limit'], 7,color=mycolor)
                ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
                ax.annotate('VCSS', xy=(.31,.12), xytext=(.33,0.12), color=mycolor, 
                    fontsize=12, fontweight='light',xycoords='axes fraction' )
        '''

    if axr <4:
        ax.set_xticklabels([])
        ax.set_xlabel('')
    
    if axc >0:
        ax.set_ylabel('')
    
    if axind ==nplots-1 or ii == len(final_unq_jvla)-1:
        plt.subplots_adjust(wspace=0.1, hspace=0.1, top = .98, bottom=0.05, 
                            left=0.05, right=0.99)
        
        if ii ==len(final_unq_jvla)-1:
            for jj in range(axr*ncols+axc+1,nrows*ncols):
                axr = int(jj/ncols)
                axc = jj - axr*ncols
                ax = axs[axr, axc]
                ax.remove()
        
        plt.savefig('SED_plots_v6_p'+str(npage)+'.pdf')
        plt.close()
        npage += 1
        
        

###############################################################
###############################################################
# Adding upper limit information in the flux table
###############################################################
###############################################################
###############################################################

tgss_lim = [np.nan]*len(wise_phot)
wenss_lim = [np.nan]*len(wise_phot)
gleam_lim = [np.nan]*len(wise_phot)
vlass_lim = [np.nan]*len(wise_phot)
vlssr_lim = [np.nan]*len(wise_phot)


for jj, row in enumerate(final_unq_jvla): 
    source = row['Source_name']
    wise_row = wise_phot.loc['Source_name',source]
    wisen = wise_row['WISEname']
    al_row = al_calcs.loc[source]
    sp_row = sp_class.loc[source]
    spshape = sp_row['Modified_Sp_Shape']
    morph = al_row['Xband_morph']
    ii = wise_phot.loc_indices[source]
    curved_cats = ['CV', 'CV:', 'PK', 'PK:']
    model_name = 'PL'
    if np.isin(spshape, curved_cats):
        model_name = 'CPL'
    
    scat = wise_row
    glmcat = atgcat.loc[wisen]
    
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
    
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
    OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,jvla_BX)
    freq_arr, flux_arr, eflux_arr, labels = modify_data_plot(source,freq_arr,
                        flux_arr, eflux_arr, labels,sp_row, wise_row )        
    
    freq_arr, flux_arr, eflux_arr, labels = VCSS_data_cleanup(source, wise_row,
                            sp_row, labels, freq_arr, flux_arr, eflux_arr)
    
    if 'VLASS' in labels:  
         vlass_lim[ii] = False
    if sp_row['VLASS_Limit']!='U':
            vlass_lim[ii] = True
    
                  
    if 'TGSS' in labels:
        tgss_lim[ii] = False
    if sp_row['TGSS_Limit'] =='U':
        tgss_lim[ii] = True
            

    if 'VLSSr' in labels:
        vlssr_lim[ii] = False
    if sp_row['VLSSr_Limit'] =='U':
        vlssr_lim[ii] = True
            

        
    if 'WENSS' in labels:
        wenss_lim[ii] = False
    if sp_row['WENSS_Limit'] =='U':   
        wenss_lim[ii] = True
            

        
    if 'GLEAM' in labels:
        gleam_lim[ii] = False

    if sp_row['GLEAM_Limit'] =='U':   
        gleam_lim[ii] = True

#
for col , cat in zip([tgss_lim, gleam_lim, vlssr_lim, wenss_lim,   vlass_lim], 
                     ['TGSS', 'GLEAM', 'VLSSr', 'WENSS', 'VLASS']):
    ind = myfdf.columns.get_loc("E"+cat)
    myfdf.insert(column='L'+cat, loc=ind+1, value=col)
        

table4 = Table.read('Tab4_Flux_measurements_spectra.csv')

  
nu_mod = np.logspace(-1.3,1.3, 100 )
nrows = 5
ncols = 3
filecounter = 0
myfilecounter = []
subplotcounter = 0
mysubplotcounter = []
plotarr = [np.nan]*(nrows*ncols)
npage = 0
nplots = nrows*ncols
plt.ion()
rgb = sns.color_palette("Set1", 10)


for ii, row in enumerate(final_unq_jvla[:1]): 
    source = row['Source_name']
    wise_row = wise_phot.loc['Source_name',source]
    wisen = wise_row['WISEname']
    al_row = al_calcs.loc[source]
    sp_row = sp_class.loc[source]
    spshape = al_row['Final_SP_Class']
    morph = al_row['Xband_morph']
    
    curved_cats = ['CV', 'CV:', 'PK', 'PK:']
    model_name = 'PL'
    if np.isin(spshape, curved_cats):
        model_name = 'CPL'
    
    
    scat = wise_row
    glmcat = atgcat.loc[wisen]
    
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
    
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
    OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,jvla_BX)
    freq_arr, flux_arr, eflux_arr, labels = modify_data_plot(source,freq_arr,
                        flux_arr, eflux_arr, labels,sp_row, wise_row )        
    
    freq_arr, flux_arr, eflux_arr, labels = VCSS_data_cleanup(source, wise_row,
                            sp_row, labels, freq_arr, flux_arr, eflux_arr)
    
    if 'VLASS' in labels and sp_row['VLASS_Limit']!='U':
        ind = labels.index('VLASS')
        vflux = flux_arr[ind]
        eflux = eflux_arr[ind]
        if eflux<0.2*vflux:
            eflux_arr[ind] = np.sqrt(eflux**2 + (vflux/5)**2)        
    
    
    name = scat['WISEname']  #Name of the object.
    msk = [freq_arr>0]
    nu_arr = freq_arr[tuple(msk)]
    s_nu= flux_arr[tuple(msk)]
    e_s_nu = eflux_arr[tuple(msk)]
    survey_faint = ['TGSS', 'WENSS', 'GLEAM', 'SUMSS', 'VLSSr']   
    freq_faint = [0.15, 0.325, 0.200, 0.843, 0.074]
    for survey, freq in zip(survey_faint, freq_faint):
        flag = sp_row[survey+'_Limit']
        if flag == 'T' or flag=='F':
            if survey not in labels:
                f_flux, f_eflux = mcsim_fit.get_faint_detection(source, survey,10)
                if f_flux>0:
                    labels.append(survey)
                    nu_arr = np.append(nu_arr, freq)
                    s_nu = np.append(s_nu, f_flux)
                    e_s_nu = np.append(e_s_nu, f_eflux)
    
    nreg = max(1, len(jvla_AX), len(jvla_BX))
    if source == 'J2325-04': nreg = 1
    
    ymin = min(flux_arr)
    ymax = max(flux_arr)
    
    axind = ii%nplots
    axr = int(axind/ncols)
    axc = axind - axr*ncols
    print(axind, axr, axc)
    #if axind ==0:
    fig,ax = plt.subplots(1,1, num=npage, figsize = (7,5))     
    #ax = axs[axr, axc]
    
    
    if 'AX' in labels and len(alpha_AX)>0:
        if np.any(jvla_AX['snr']>50) and alpha_AX[1]>-999:
            mind = labels.index('AX')
            FAX_10GHZ = flux_arr[mind]
            EAX_10GHZ = eflux_arr[mind]
            plot_bowtie(10,alpha_AX, FAX_10GHZ, EAX_10GHZ,ax, zorder=9)
    if 'BX' in labels and len(alpha_BX)>0:
        if np.any(jvla_BX['snr']>50) and alpha_BX[1]>-999:
            mind = labels.index('BX')
            FBX_10GHZ = flux_arr[mind]
            EBX_10GHZ = eflux_arr[mind]
            plot_bowtie(10,alpha_BX, FBX_10GHZ, EBX_10GHZ,ax, zorder=9)           
    
    if sp_row['TGSS_Limit'] =='U':
        noise = wise_row['TGSS_Limit']
        plot_upper_limit(0.15, 2*noise, noise, 9)
        ymin = min(ymin, noise)
        ymax = max(ymax, 3*noise)
    
    if sp_row['VLSSr_Limit'] =='U':
        labels.append('VLSSr')
        noise = wise_row['VLSSr_Limit']
        plot_upper_limit(0.074, 2*noise, noise, 7)
        ymin = min(ymin, noise)
        ymax = max(ymax, 3*noise)
    
    if sp_row['WENSS_Limit'] =='U':            
        noise = wise_row['WENSS_Limit']
        if noise>0:
            labels.append('WENSS')
            plot_upper_limit(0.325, 2*noise, noise, 7)
            ymin = min(ymin, noise)
            ymax = max(ymax, 3*noise)
    
    if sp_row['GLEAM_Limit'] =='U':            
        noise = wise_row['GLEAM_Limit']
        if noise>0:
            labels.append('GLEAM')
            plot_upper_limit(0.200, 2*noise, noise, 7)
            ymin = min(ymin, noise)
            ymax = max(ymax, 3*noise)
    
    ax.errorbar(nu_arr, s_nu, yerr=e_s_nu, marker ='o', 
                linestyle='none', capsize=3, markersize = 4, color='k', 
                elinewidth = 1, zorder = 10)            
    
    
    ax.set_title(source,  fontsize=16, fontweight='medium')
    ax.set_xscale('log')
    ax.set_yscale('log')
    if spshape=='PL':
        ax.set_ylim(ymin/2, 3*ymax)
    else:
        ax.set_ylim(ymin/2, 3*ymax)
    ax.set_xlim(0.03, 20)
    ax.set_xlabel(r'log ($\nu / GHz)$', fontdict=dict(fontsize=14, ))
    ax.set_ylabel(r'log ($S_{\nu}/mJy)$', fontdict=dict(fontsize=14))
    ax.tick_params(axis = 'both', which = 'minor', direction='in', 
                   length=4, top=True, right=True)
    ax.tick_params(axis = 'both', which = 'major', direction='in',
                   length=9, top=True, right=True, labelsize=14)
    
    
    #ax.text(0.9,0.85, spshape, 
    #        transform=ax.transAxes, color='k', 
    #        fontsize = 16,bbox = dict(boxstyle='round',fc='w', ec='k'))
    #ax.text(0.9,0.1, morph, 
    #        transform=ax.transAxes, color='k', 
    #        fontsize = 16,bbox = dict(boxstyle='round',fc='w', ec='k'))
    
    
    fitRes = fit_results_ext_tab.loc[source]
    par_pl = fitRes['PL_s0', 'PL_al']
    par_cpl = fitRes['CPL_s0', 'CPL_al', 'CPL_nup']
    
    nu_mod = np.logspace(-3, 3.0, 5000)
    
    s1 = rmfit.power_law(nu_mod, *par_pl)
    s2 = rmfit.CPL(nu_mod, *par_cpl)
    lincols = ['#66c2a5', '#fc8d62','#8da0cb','#e78ac3','#a6d854']
    
    if spshape in ['PL', 'PL:', 'I', 'F']:
        if spshape in ['PL', 'PL:']:
            l1, = ax.plot(nu_mod, s1, linestyle = 'dotted', label='PL', color=lincols[0], 
                  lw=1, zorder = 1)       
    else:
        l2, = ax.plot(nu_mod, s2, linestyle = 'dotted', label='CV', color=lincols[1], 
                  lw=1, zorder = 1)
    
    z = al_row['redshift']
    if z>0:
        ax.text(0.03, 0.9,'z={0}'.format(z), fontsize=14, bbox=dict(facecolor='white', alpha=0.5, 
                                  boxstyle='round', ec = 'gray'), 
                        transform=ax.transAxes, zorder = 11)
   
    
    ax.axvspan(1.37, 1.43,facecolor=rgb[0], alpha=0.2, zorder=1 )
    ax.axvspan(8, 12,facecolor=rgb[3], alpha=0.2, zorder=1 )
    ax.axvspan(.139, .156,facecolor=rgb[4], alpha=0.2, zorder=1 )
    ax.axvspan(2, 4,facecolor=rgb[8], alpha=0.2, zorder=1 )
    

    ax.annotate('NVSS', xy=(.51,.2), xytext=(.55,0.2), color=rgb[0], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    ax.annotate('TGSS', xy=(.13,.2), xytext=(.20,0.2), color=rgb[4], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )    
    ax.annotate('VLA-X', xy=(.83,.2), xytext=(.85,.2), color=rgb[3], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    if (wisen =='040403.61-243600.1') or (wisen == '061200.23-062209.1'):
        ax.annotate('VLASS', xy=(.64,.1), xytext=(.67,.1), color=rgb[8], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    else:
        ax.annotate('VLASS', xy=(.64,.12), xytext=(.65,.12), color=rgb[8], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    
    if 'VLITE' in labels:
        mycolor = 'tab:red'
        #ax.errorbar(vlite_freq, vlite_flux, yerr=vlite_eflux, marker ='*', 
        #    linestyle='none', capsize=3, markersize = 4, color=mycolor, 
        #    elinewidth = 1, zorder = 10)

        ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
        ax.annotate('VCSS', xy=(.31,.2), xytext=(.33,0.2), color=mycolor, 
            fontsize=12, fontweight='light',xycoords='axes fraction' )   
    
    if ii not in [103,139]:
        if 'WENSS' in labels or 'TEXAS' in labels or 'WISH' in labels:
            ax.axvspan(.322, .328,facecolor=rgb[2], alpha=0.2 )
        #ax.annotate('VLITE', xy=(.28,.2), xytext=(.28,0.2), color=rgb[2], 
        #        fontsize=12, fontweight='medium',xycoords='axes fraction' )
    else:
        ax.axvspan(.35, .38,facecolor=rgb[2], alpha=0.2, zorder=1 )
   
    if 'GLEAM' in labels: 
        ax.axvspan(0.17, 0.23,facecolor=rgb[6], alpha=0.1, zorder=1 )
        ax.annotate('GLEAM', xy=(.16,.12), xytext=(.23,0.12), color=rgb[6], 
                fontsize=12, fontweight='light',xycoords='axes fraction' ) 
    if 'WENSS' in labels:
        ax.annotate('WENSS', xy=(.31,.12), xytext=(.31,0.12), color=rgb[2], 
                fontsize=12, fontweight='light',xycoords='axes fraction' ) 
    if 'VLSSr' in labels:
        ax.axvspan(.073, .075,facecolor=rgb[1], alpha=0.1, zorder=1)
        ax.annotate('VLSSr', xy=(.05,.12), xytext=(.06,0.12), color=rgb[1], 
                fontsize=12, fontweight='light',xycoords='axes fraction' )
    if source =='J1651+34':
        ax.axvspan(3.85, 5.85,facecolor='gray', alpha=0.1, zorder=1 )
        ax.annotate('GB6', xy=(.75,.2), xytext=(.75,0.2), color='gray', 
                fontsize=12, fontweight='light',xycoords='axes fraction' )
    if 'RACS' in labels:
        ax.axvspan(0.743, 1.031,facecolor=rgb[8], alpha=0.1, zorder=1 )
        ax.annotate('RACS', xy=(.46,.12), xytext=(.46,0.12), color=rgb[8], 
                    fontsize=12, fontweight='medium',xycoords='axes fraction' )
       
    if 'Texas' in labels:
        ax.annotate('TEXAS', xy=(.32,.2), xytext=(.33,0.2), color=rgb[2], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    '''
        else:    
            if sp_row['VLITE_Limit']=='F':               
                mycolor= 'tab:cyan'
                ax.errorbar(vlite_freq, vlite_flux, yerr=vlite_eflux, marker ='*', 
                    linestyle='none', capsize=3, markersize = 4, color=mycolor, 
                    elinewidth = 1, zorder = 10)
                ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
                ax.annotate('VCSS', xy=(.31,.12), xytext=(.33,0.12), color=mycolor, 
                    fontsize=12, fontweight='light',xycoords='axes fraction' )

             
            elif wise_row['VLITE_Limit']>0:
                mycolor = 'tab:cyan'
                plot_upper_limit(0.338, 2.5*wise_row['VLITE_Limit'], 
                        1.5*wise_row['VLITE_Limit'], 7,color=mycolor)
                ax.axvspan(0.32, 0.368,facecolor=mycolor, alpha=0.1, zorder=1 )
                ax.annotate('VCSS', xy=(.31,.12), xytext=(.33,0.12), color=mycolor, 
                    fontsize=12, fontweight='light',xycoords='axes fraction' )
        

    if axr <4:
        ax.set_xticklabels([])
        ax.set_xlabel('')
    
    if axc >0:
        ax.set_ylabel('')
    
    #if axind ==nplots-1 or ii == len(final_unq_jvla)-1:
    #    plt.subplots_adjust(wspace=0.1, hspace=0.1, top = .98, bottom=0.05, 
                            left=0.05, right=0.99)
        
    #    if ii ==len(final_unq_jvla)-1:
    #        for jj in range(axr*ncols+axc+1,nrows*ncols):
    #            axr = int(jj/ncols)
    #            axc = jj - axr*ncols
    #            ax = axs[axr, axc]
    #            ax.remove()
    '''    
    plt.tight_layout()
    plt.savefig('RadioSED_plots/'+source+'_SED_plots_v3_p'+str(npage)+'.pdf')
    plt.close()
    npage += 1
        


# Frequency vs sensitivity plot

df = pd.DataFrame({'survey' : ['VLA-A/B', '', 'VLASS',  'NVSS', 'FIRST', 'SUMSS', 'WENSS', 
          'GLEAM', 'TGSS', 'VLSSr', 'VCSS', 'RACS'], 
'cen_freq': [10, 10, 3, 1.4, 1.4, 0.843, 0.325, 0.200, 0.150, 0.074, 0.338, 0.887],
'sensitivity' : [0.013, 0.013, 0.15, 0.45, 0.15, 1.1, 3., 8, 3.5, 100, 5, 0.3],
'res': [0.2, 0.6, 2.5, 45, 5, 45, 54, 100, 25, 75, 10, 15], 
'bw': [2, 2, 1, 0.025, 0.025,0.002, 0.003, 0.030,0.0085,  0.001, 0.030,
       0.144], 
 'txt_off':[(-20,5), (-20,5), (-20,8), (-20,15), (-20,-18), (-20,18), 
            (-20,-30),(-20, 25), (-20,-23), 
            (-20,-33), (-20,8), (-20,-20)]})

labeldict = {'size':14, 'weight':'normal', 'stretch':'normal', 'family':'serif', 
             'style':'normal', 'variant':'normal'}
fig = plt.figure()
ax = sns.scatterplot(data = df, x = df.cen_freq, y = df.sensitivity, hue=df.res, size=df.res, sizes=(10, 1000), 
                     alpha=0.5, edgecolor='k')
for line in range(0,df.shape[0]):
     ax.annotate(df.survey[line], xy = (df.cen_freq[line], df.sensitivity[line]),
                 xytext= df.txt_off[line], textcoords='offset pixels', 
             size=10, color='k', weight='medium', zorder=4, alpha=0.9)
     #ax.text(df.cen_freq[line], df.sensitivity[line], df.survey[line], horizontalalignment='center', 
     #        size=10, color='k', weight='medium', zorder=4, alpha=0.9)
        
for line in range(0,df.shape[0]):
    ax.errorbar(df.cen_freq[line], df.sensitivity[line], xerr= df.bw[line],marker ='.', 
                    linestyle='none', markersize = 1, capsize=3, color='grey', 
                    elinewidth = 1, )
ax.legend_.remove()
ax.set_ylim(0.005, 200)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.05, 20)
ax.minorticks_on()
ax.tick_params(axis= 'both', which='major', labelsize=12, length = 7, direction = 'in',top = True, 
                right = True)
ax.tick_params(axis= 'both', which='minor', length = 4, direction = 'in',top = True, 
                right = True)

ax.set_ylabel('Sensitivity (mJy/beam)', **labeldict)
ax.set_xlabel('Frequency (GHz)',  **labeldict)

#add legend
from matplotlib.patches import Patch
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

mycmap = plt.cm.plasma
norm = Normalize(vmin= min(df.res), vmax = max(df.res))
mycolor = mycmap(norm(df.res[0]))
patch1 = Line2D([0], [0], marker='o', color='grey', label='0.5"', alpha= 0.5, 
                          markerfacecolor=mycolor, markersize=3, linestyle='none')
mycolor = mycmap(norm(df.res[4]))
patch2 = Line2D([0], [0], marker='o', color='grey', label='5"',alpha= 0.5, 
                          markerfacecolor=mycolor, markersize=9, linestyle='none')
mycolor = mycmap(norm(df.res[3]))
patch3 = Line2D([0], [0], marker='o', color='grey', label='50"',alpha= 0.5, 
                          markerfacecolor=mycolor, markersize=20, linestyle='none')

leg = plt.legend(handles=[patch1, patch2, patch3], loc='upper right', 
                  handlelength=4, labelspacing = 1)
 
ax.tick_params(axis= 'both', which='major', labelsize=14, length = 7, labeltop = False, 
                    labelbottom = True, direction = 'in', right  = True, top=True)
ax.tick_params(axis= 'both', which='minor', labelsize=14, length = 4, labeltop = False, 
                    labelbottom = True, direction = 'in', right = True, top=True)
ax.minorticks_on()
 
plt.tight_layout()

plt.savefig('surveys_labels_v1.pdf', dpi=500)



'''
###############################################################################
##############Fitiing for sources with Peaked or Curved SED:
###############################################################################
###############################################################################

'''        
'''
sources =  []
spcodes = ['CV', 'CV:', 'PK', 'PK:']
for row in sp_class:
    if row['Clean_Code']=='T' and row['Final_SP_Class'] in spcodes:
        sources.append(row['Source_name'])
                          
        
        
models = ['GCV', 'EFFA', 'IFFA', 'SSA']
colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
#rgb = sns.xkcd_palette(colors)
rgb = sns.color_palette("Set3", 6)
sns.set(style="ticks", palette="pastel", color_codes=True)
plt.rcParams.update({'lines.markeredgewidth': 1})
lincols = ['#66c2a5', '#fc8d62','#8da0cb','#e78ac3','#a6d854']

plot_conf = False
plot_parreg = False

#for i in range(0,155): 
for ii, source in enumerate(sources):
    wise_row = wise_phot.loc[source]
    wisen = wise_row['WISEname']
    sp_row = sp_class.loc[source]
    morph = sp_row['Xband_morph']
    lomorph = sp_row['Visual_morph']
    spshape = sp_row['Final_SP_Class']
    mc_par = par_fits[source]
    fig, ax= plt.subplots(1,1,figsize= (6.52, 8.34))
    scat = wise_row
    glmcat = atgcat.loc[wisen]
    
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
         
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
        OIR, sp_flux, labels = fprep.data_prep(wisen, glmcat, scat, jvla_AX,jvla_BX)
    freq_arr, flux_arr, eflux_arr, labels = modify_data_plot(source,freq_arr,
                            flux_arr, eflux_arr, labels,sp_row, wise_row )        
        
        
        
    name = scat['WISEname']  #Name of the object.
    msk = [freq_arr>0]
    nu_arr = freq_arr[tuple(msk)]
    s_nu= flux_arr[tuple(msk)]
    e_s_nu = eflux_arr[tuple(msk)]
    if 'AX' in labels and len(alpha_AX)>0:
        if np.any(jvla_AX['snr']>50) and alpha_AX[1]>-999:
            mind = labels.index('AX')
            FAX_10GHZ = flux_arr[mind]
            EAX_10GHZ = eflux_arr[mind]
            plot_bowtie(10,alpha_AX, FAX_10GHZ, EAX_10GHZ,ax, zorder=9)
    if 'BX' in labels and len(alpha_BX)>0:
        if np.any(jvla_BX['snr']>50) and alpha_BX[1]>-999:
            mind = labels.index('BX')
            FBX_10GHZ = flux_arr[mind]
            EBX_10GHZ = eflux_arr[mind]
            plot_bowtie(10,alpha_BX, FBX_10GHZ, EBX_10GHZ,ax, zorder=9)
            
                
    if sp_row['TGSS_Limit'] =='U':
        noise = wise_row['TGSS_Limit']
        plot_upper_limit(0.15, 2*noise, noise, 9)
        
    if sp_row['VLSSr_Limit'] =='U':
        labels.append('VLSSr')
        noise = wise_row['VLSSr_Limit']
        plot_upper_limit(0.074, 2*noise, noise, 7)

    if sp_row['WENSS_Limit'] =='U':            
        noise = wise_row['WENSS_Limit']
        if noise>0:
            labels.append('WENSS')
            plot_upper_limit(0.325, 2*noise, noise, 7)
            
    if sp_row['GLEAM_Limit'] =='U':            
        noise = wise_row['GLEAM_Limit']
        if noise>0:
            labels.append('GLEAM')
            plot_upper_limit(0.200, 2*noise, noise, 7)
    survey_faint = ['TGSS', 'WENSS', 'GLEAM', 'SUMSS', 'VLSSr']   
    freq_faint = [0.15, 0.325, 0.200, 0.843, 0.074]
    for survey, freq in zip(survey_faint, freq_faint):
        flag = sp_row[survey+'_Limit']
        if flag == 'T' or flag=='F':
            if survey not in labels:
                f_flux, f_eflux = mcsim_fit.get_faint_detection(source, survey,10)
                if f_flux>0:
                    labels.append(survey)
                    nu_arr = np.append(nu_arr, freq)
                    s_nu = np.append(s_nu, f_flux)
                    e_s_nu = np.append(e_s_nu, f_eflux)
                
    ax.errorbar(nu_arr, s_nu, yerr=e_s_nu, marker ='o', 
                linestyle='none', capsize=3, markersize = 4, color='k', 
                elinewidth = 1, zorder = 10)
    nreg = max(1, len(jvla_AX), len(jvla_BX))
    if source == 'J2325-04': nreg = 1
    # Adding individual points for multi-coomponent sources.
    config = sp_row['Final_array']
    mystr = '{0}\nz: {1}\nX Morph: {2}\nLo Morph: {3}\nSED: {4}\nconfig: {5}\n'.format(
                           source,sp_row['redshift'], morph, lomorph, spshape, config )
    mystr = spshape
    if sp_row['redshift'] >0:
        mystr = mystr+', z: {0}'.format(sp_row['redshift'])
    ax.text(0.03, 0.9,mystr, fontsize=14, bbox=dict(facecolor='white', alpha=0.5, 
                                      boxstyle='round', ec = 'gray'), 
                            transform=ax.transAxes, zorder = 11)
            
           
    # Adjusting the labels, axes scales and range.  
    ax.set_title(source,  fontsize=16, fontweight='medium')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(0.01*min(flux_arr), 10*max(flux_arr))
    ax.set_xlim(0.03, 20)
    ax.set_xlabel(r'log ($\nu / GHz)$', fontdict=dict(fontsize=14, ))
    ax.set_ylabel(r'log ($S_{\nu}/mJy)$', fontdict=dict(fontsize=14))
    ax.tick_params(axis = 'both', which = 'minor', direction='in', 
                   length=4, top=True, right=True)
    ax.tick_params(axis = 'both', which = 'major', direction='in',
                   length=9, top=True, right=True, labelsize=14)

        
            
    nu_mod = np.logspace(-3, 3.0, 5000)
    
    fitRes = fit_results_tab.loc[source]
    par_effa = fitRes['EFFA_s0', 'EFFA_al', 'EFFA_nup']
    par_iffa = fitRes['IFFA_s0', 'IFFA_al', 'IFFA_nup']
    par_ssa = fitRes['SSA_s0', 'SSA_al', 'SSA_nup']
                
                
    s3 = rmfit.EFFA_func(nu_mod, *par_effa )
    s4 = rmfit.IFFA_func(nu_mod, *par_iffa)
    s5 = rmfit.SSA_func(nu_mod, *par_ssa)
        
    pars, gc_fit_val, tcarlo  = mcsim_fit.mcsim_fitting(wise_row, 'GCV', atgcat, vla_ax_grp, 
                                                            vla_bx_grp, 1000)
    gc_fit =  dict.fromkeys(['sp', 'al_thin', 'al_thick', 'nupeak' ])
    for key, prow in zip(cpl_fit.keys(), cpl_fit_val):
        pardict = {} 
        for k, v in zip(pars, prow):
            pardict[k] = v
        cpl_fit[key] = pardict            
        
    for p1, p2, p3, p4 in zip(tcarlo[0][:],tcarlo[1][:],tcarlo[2][:] ,tcarlo[3][:]):
        s7 = rmfit.gen_curve_func(nu_mod, *[p1, p2, p3, p4])   
        ax.plot(nu_mod, s7, linestyle = '-', color='gray', 
                  lw=0.5, alpha = 0.1, zorder = 1)
    med_pars = [np.median(tcarlo[0][:]),np.median(tcarlo[1][:]),np.median(tcarlo[2][:]), 
                     np.median(tcarlo[3][:])]
    med_epar = [np.std(tcarlo[0][:]),np.std(tcarlo[1][:]),np.std(tcarlo[2][:]), 
                np.median(tcarlo[3][:])]
        
        
    s6 = rmfit.gen_curve_func(nu_mod, *med_pars) 
    l6, = ax.plot(nu_mod, s6, linestyle = '-', label='GCV-MC', color=lincols[1], 
                      lw=2, zorder = 3)
        
    l3, = ax.plot(nu_mod, s3, linestyle = '-.', label='EFFA', color=lincols[2], 
                  lw=2, zorder = 3)
    l4, = ax.plot(nu_mod, s4, linestyle = 'dotted', label='IFFA',  color=lincols[3], 
                  lw=2, zorder= 3)
    l5, = ax.plot(nu_mod, s5, linestyle = ':', label='SSA',  color=lincols[4], 
                  lw=2, zorder = 3)
        
    cpl_al = par_cpl['CPL_al']
    cpl_q = par_cpl['CPL_nup']
    cpl_nup = np.exp(-cpl_al/(2*cpl_q))
    cpl_nup_err = abs((-cpl_nup/(2*cpl_q))*(perr_cpl[1] + (cpl_al*perr_cpl[2]/cpl_q)))

    if plot_parreg:
        if abs(med_pars[2])>0.1:
            lpert = 15.9
            hpert = 84.1
            mcnup = np.exp(-tcarlo[1][:]/(2*tcarlo[2][:]))
            mcsp = rmfit.CPL(mcnup, tcarlo[0][:], tcarlo[1][:], tcarlo[2][:])
            mc_low = [np.percentile(mcsp, lpert),
                           np.percentile(mcnup, lpert)]
            mc_hig = [np.percentile(mcsp, hpert),
                           np.percentile(mcnup, hpert)]
            
            x = np.median(mcsp)-mc_low[0]
            y = np.median(mcnup)-mc_low[1]
            dx = mc_low[0]+mc_hig[0]
            dy = mc_low[1]+mc_hig[1]
            
            #x = cpl_fit['speak']['Median'] - cpl_fit['speak']['e1Low']
            #x = cpl_fit['speak']['Median'] - cpl_fit['speak']['e1Low']

            rect_cpl = Rectangle((y,x),dy, dx, fc = lincols[1], 
                               label='CV-MC', alpha = 0.5,ec='white', linewidth=2.0 )            
            try:
                ax.add_patch(rect_cpl)
            except: continue
        
    #par_str = r' $\alpha_{PL}:$'+' {0:.2f},'.format(par_pl['PL_al']) 
    #par_str = par_str + r' $\alpha_{CV}:$'+' {0:.2f}'.format(par_cpl['CPL_al']) +r'$\pm$'+' {0:.2f}\n'.format(perr_cpl[1])
    #par_str = par_str+r' $\nu^{p}_{CV}:$'+' {0:.1e}'.format(cpl_nup)+r'$\pm$'+' {0:.1e} GHz\n'.format(cpl_nup_err)
    #par_str = par_str+r' $q_{CV}:$'+' {0:.1f}'.format(cpl_q)+r'$\pm$'+' {0:.1f}'.format(perr_cpl[2])
    #par_str = par_str+r' $S_{0,CV}:$'+' {0:.1f}'.format(par_cpl[0])+r'$\pm$'+' {0:.1f}\n'.format(perr_cpl[0])
    #par_str = par_str + r' $\alpha_{CV,M}:$'+' {0:.2f}'.format(med_pars[1])+r'$\pm$'+' {0:.2f},'.format(med_epar[1])
    #par_str = par_str+r' $q_{CV,M}:$'+' {0:.1f}\n'.format(med_pars[2])+r'$\pm$'+' {0:.2f},\n'.format(med_epar[2])
    #par_str = par_str+r' $S_{0,CV,M}:$'+' {0:.1f}'.format(med_pars[0])+r'$\pm$'+' {0:.2f}'.format(med_epar[0])
    

    ax.text(0.27, 0.80,par_str, fontsize=8, bbox=dict(facecolor='white', alpha=0.5, 
                                      boxstyle='round', ec = 'gray'), 
                            transform=ax.transAxes, zorder=12)
    
   
    ax.axvspan(1, 2,facecolor=rgb[0], alpha=0.2, zorder=1 )
    ax.axvspan(8, 12,facecolor=rgb[3], alpha=0.2, zorder=1 )
    ax.axvspan(.12, .18,facecolor=rgb[4], alpha=0.2, zorder=1 )
    ax.axvspan(2.5, 4,facecolor=rgb[5], alpha=0.2, zorder=1 )
    

    ax.annotate('NVSS', xy=(.51,.2), xytext=(.55,0.2), color=rgb[0], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    ax.annotate('TGSS', xy=(.13,.2), xytext=(.20,0.2), color=rgb[4], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )    
    ax.annotate('VLA-X', xy=(.83,.2), xytext=(.85,.2), color=rgb[3], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    if (wisen =='040403.61-243600.1') or (wisen == '061200.23-062209.1'):
        ax.annotate('VLASS', xy=(.64,.1), xytext=(.67,.1), color=rgb[5], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    else:
        ax.annotate('VLASS', xy=(.64,.2), xytext=(.67,.2), color=rgb[5], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )
    
    if ii not in [103,139]:
        ax.axvspan(.296, .4,facecolor=rgb[2], alpha=0.2 )
        #ax.annotate('VLITE', xy=(.28,.2), xytext=(.28,0.2), color=rgb[2], 
        #        fontsize=12, fontweight='medium',xycoords='axes fraction' )
    else:
        ax.axvspan(.35, .38,facecolor=rgb[2], alpha=0.2, zorder=1 )
   
    if 'GLEAM' in labels: 
        ax.annotate('GLEAM', xy=(.16,.12), xytext=(.23,0.12), color=rgb[4], 
                fontsize=12, fontweight='light',xycoords='axes fraction' ) 
    if 'WENSS' in labels:
        ax.annotate('WENSS', xy=(.31,.12), xytext=(.33,0.12), color=rgb[2], 
                fontsize=12, fontweight='light',xycoords='axes fraction' ) 
    if 'VLSSr' in labels:
        ax.axvspan(.073, .075,facecolor=rgb[1], alpha=0.1, zorder=1)
        ax.annotate('VLSSr', xy=(.02,.2), xytext=(.02,0.2), color=rgb[4], 
                fontsize=12, fontweight='light',xycoords='axes fraction' )
    if source =='J1651+34':
        ax.axvspan(3.85, 5.85,facecolor='gray', alpha=0.1, zorder=1 )
        ax.annotate('GB6', xy=(.75,.2), xytext=(.75,0.2), color='gray', 
                fontsize=12, fontweight='light',xycoords='axes fraction' )

    if 'Texas' in labels:
        ax.annotate('TEXAS', xy=(.32,.2), xytext=(.33,0.2), color=rgb[2], 
                fontsize=12, fontweight='medium',xycoords='axes fraction' )

    
    #ax.annotate(name,xy=(.3,1.1), xytext=(.3,1.1), color='k', 
    #         fontsize=12, fontweight='light' ,xycoords='axes fraction' )
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.subplots_adjust(bottom = 0.45, right = 0.99)
    ax3 = plt.axes([0.02, 0.01, 0.97, 0.35])
    ax3.set_xticks([])
    ax3.set_yticks([])
    imgfile = 'figures/'+ source.strip()+'_comb_v3.png'
    if os.path.exists(imgfile):
        img = mpimg.imread(imgfile)
        ax3.imshow(img, aspect='auto')
        ax3.set_xticks([])
        ax3.set_yticks([])
    plt.savefig('RadioSED_plots/peak_sources/'+name+'_rsed_9.pdf', dpi = 300)
    plt.close()
        
   
#Combine Before and after
import glob 
myimgs = glob.glob('./RadioSED_plots/New_Fit_Color_Plots/*_rsed_11.pdf')
import matplotlib.image as mpimg
for im in myimgs[:1]:
    base = os.path.basename(im)
    prefix = base.split('_')[0]
    prim_image = 'RadioSED_plots/New_Fit_Color_Plots/'+prefix+'_rsed_9.pdf'
    secn_image = 'RadioSED_plots/New_Fit_Color_Plots/'+prefix+'_rsed_11.pdf'
    fig, ax= plt.subplots(1,2,figsize= (10,10))
    ax1 = ax[0] 
    ax2 = ax[1]
    if os.path.exists(prim_image):
        img = mpimg.imread(prim_image)
        ax1.imshow(img, aspect='auto')
        ax1.set_xticks([])
        ax1.set_yticks([])
    if os.path.exists(secn_image):
        img = mpimg.imread(secn_image)
        ax2.imshow(img, aspect='auto')
        ax2.set_xticks([])
        ax2.set_yticks([])
    plt.tight_layout() 
    plt.savefig('RadioSED_plot/'+source+'_comb_v4.png', dpi = 500)
    plt.close()
'''