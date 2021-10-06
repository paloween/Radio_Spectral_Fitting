#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 23:45:58 2019

@author: pallavipatil
"""

#Interactive SED fitting

import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor
from astropy.io import fits,ascii
import matplotlib.pyplot as plt
from astropy.table import Table, unique
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
import aplpy as apl
from astropy.visualization import ZScaleInterval
import scipy
import pandas as pd
from Radio_Models_func import RadioModelsFits
from file_prep_radio_fitting import LoadTables, FilePrepUtils



##Accessing all different data files and creating astropy tables.

loadt = LoadTables()
fprep = FilePrepUtils()
rmfit = RadioModelsFits()
atscat, atgcat, vla_ax_grp, vla_bx_grp = loadt.get_tabs()
sp_class = loadt.get_sp_info()



plots_dir = './Plots_all/'
fitsdir_AX = '/Users/pallavipatil/Desktop/VLA/VLA-NewWork/AX/fitsfiles/'
fitsdir_BX = '/Users/pallavipatil/Desktop/VLA/VLA-NewWork/BX/fitsfiles/'

scan_info = pd.ExcelFile('./Datasets/JVLA_Scans_info.xlsx')
dfsinfo = pd.read_excel(scan_info, 'Summary')
scan_info.close()
ax_list = pd.Series(dfsinfo['12B-217'])
bx_list = pd.Series(dfsinfo['12A-064'])



'''##########################################################################'''

'''
This part sets up the figure area, boxes, subplots and buttons.
'''
'''##########################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
'''

#fig = plt.figure(figsize = (12,6))
#grid = plt.GridSpec(1,4, wspace = 0.2)
#imax = fig.add_subplot(grid[:, -1])
#ax = fig.add_subplot(grid[:,:-1])
fig, ax= plt.subplots(figsize= (12,8))
plt.subplots_adjust(left=0.20, bottom=0.23, right=0.7,top=0.8)
cursor = Cursor(ax, useblit=True, color='k', linewidth=1 )

infoax = plt.axes([0.025,0.25,0.12, 0.5], facecolor= '#C8C9C7', alpha=0.1)
infoax.tick_params(axis='both', bottom='off', left='off', labelbottom='off', labelleft='off')

# Red Chi square information:
axchi = plt.axes([0.025, 0.07, 0.12, 0.15],facecolor= '#C8C9C7', alpha=0.1)
axchi.tick_params(axis='both', bottom='off', left='off', labelbottom='off', labelleft='off')

#####Sliders
axcolor = 'lightblue' # These are slides for parametes. First, setting up an area for each slider.
axalpha = plt.axes([0.20, 0.07, 0.47, 0.02], facecolor=axcolor)
axs0 = plt.axes([0.20, 0.1, 0.47, 0.02], facecolor=axcolor)
axnu = plt.axes([0.20, 0.13, 0.47, 0.02], facecolor=axcolor)
'''
# Using a predefined Slider function from matplotlib.Widgets. You can specify 
the range for the parameter and an intial value for each slider.'''
salpha = Slider(axalpha, r'$\alpha$', -3, 3.0, valinit=-1) #Spectral index
snu = Slider(axnu, r'log $\nu_{p}/q$', -3, 2, valinit=-1)# Turnover freq
ss0 = Slider(axs0, r'log $S_{0}$', -1, 3, valinit=0) # Normalization


# Model selection
''' A radio buttion allows to select between different options.'''
rax = plt.axes([0.025, 0.80, 0.12, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ( 'EFFA', 'IFFA', 'SSA', 'PL','CPL'), active=0)

# Go to Text Box
axbox = plt.axes([0.21, 0.015, 0.07, 0.03])
fastfor= Button(axbox, 'Fward 10', color=axcolor, hovercolor = '0.975')

axbck = plt.axes([0.13, 0.015, 0.07, 0.03])
fastbck= Button(axbck, 'Bward 10', color=axcolor, hovercolor = '0.975')

#save button
'''To save the best fit parameters in a file for the analysis.''' 
axsave = plt.axes([0.29, 0.015, 0.07, 0.03])
bsave = Button(axsave, 'Save', color=axcolor, hovercolor = '0.975')

#Next and Previous button
'''Switch between different targets.'''
axprev = plt.axes([0.37, 0.015, 0.07, 0.03])
axnext = plt.axes([0.45, 0.015, 0.07, 0.03])
bnext = Button(axnext, 'Next',color=axcolor, hovercolor='0.975')
bprev = Button(axprev, 'Previous',color=axcolor, hovercolor='0.975')

# Refit button
'''If you change the guess parameters, a new fit will be generated using new 
guess values.'''
refitax = plt.axes([0.53,0.015,0.07,0.03])
ref_button = Button(refitax, 'Fit', color=axcolor, hovercolor='0.975')

#replot button
'''If you want to play around with the different parameters and see how each 
model looks like.'''
axrepl= plt.axes([0.61, 0.015, 0.07, 0.03])
brepl = Button(axrepl, 'Replot', color=axcolor, hovercolor = '0.975')

# Reset Button
# Using a predefined buttion function. 
resetax = plt.axes([0.69, 0.015, 0.07, 0.03])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


chi_ax = plt.axes([0.77,0.07,0.21,0.21])
alpha_ax = plt.axes([0.2,0.83,0.50,0.13])
alpha_ax.tick_params(axis = 'both', which = 'minor', direction='in', 
                     length=4, top=True, right=True)
alpha_ax.tick_params(axis = 'both', which = 'major', direction='in', 
                     length=9, top=True, right=True)
alpha_ax.set_xscale('log')
alpha_ax.set_xlim(0.001, 1000)
alpha_ax.set_ylim(-2, 2)
alpha_ax.axhline(0, xmin = 0.001, xmax=100,linestyle='-.')
alpha_ax.set_ylabel(r'$\alpha$', fontdict=dict(fontsize=10))  


'''##########################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
'''
   



        

#alpha_ax = plt.axes([0.2,0.83,0.50,0.13])
#chi_ax = plt.axes([0.75,0.07,0.23,0.23])

'''This file will store best fit guess parameter for each model.''' 
'''
fitTab = ascii.read('GuessPar_Radiofits.csv', format='basic', delimiter=',')
fitTab.add_index('Name')

fitRes = Table(names = ('Source_name', 'RA', 'Dec','flux_AX', 'flux_err_AX', 'SpIdx_AX', 'SpIdx_err_AX', 
                        'flux_BX', 'flux_err_BX', 'SpIdx_BX', 'SpIdx_err_BX','redshift', 'PL_S0', 
                        'PL_alpha', 'PL_rchi', 'CPL_S0', 'CPL_alpha', 'CPL_q','CPL_rchi', 'No. Points'),
                dtype = ('S20', 'S20','S20', 'f8', 'f8','f8','f8','f8','f8','f8','f8','f8','f8',
                         'f8','f8','f8','f8','f8','f8','f8'))
'''
fitRes = fit_results_tab
fitRes.add_index('source')
plot_flag = 'all'
fig_flag = True
'''A class to do everything.'''
models = ['PL', 'CPL', 'EFFA', 'IFFA', 'SSA']

class Plot_class:
    def __init__(self):
        '''Initiallizing all variables that will be used in the class functions.'''
        self.l1 = None 

        self.nu_mod = np.logspace(-3, 3.0, 5000)
        self.s1 = None
        self.fit_res = None
        self.nu_arr=None
        self.s_nu=None
        self.e_s_nu=None
        self.keypress = None
        self.clickx_data = None
        self.clicky_data = None
        self.name = None
        self.tmp=None
        self.ind = 0
        self.to = None
        self.perr_res  = None
        self.fit_pl = None
        self.fit_cpl = None
        self.perr_pl = None
        self.perr_cpl = None 


    def plotting(self,ii):
        '''This is the most important function. It plots the observations and 
        runs the non-linear chi-square miminization routine on three models.'''
        fig_flag = True
        
        self.name = atgcat['WISEname'][ii]  #Name of the object. 
        scat = atscat.loc[self.name]
        glmcat = atgcat.loc[self.name]
        sname = scat['Source_name']
        jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == self.name]
        jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == self.name]
        freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
        OIR, sp_flux, labels = fprep.data_prep(self.name, glmcat, scat, jvla_AX,
                                               jvla_BX)
        for alpha in [alpha_AX, alpha_BX, alpha_GL]:
            freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, flux_arr,
                                                               eflux_arr, alpha)
        

        ''' If multiple components are present, it will add up the flux densities
        and errors to give total value.'''
        
        FAX_10GHZ = sp_flux[0]
        EAX_10GHZ = sp_flux[1]

        FBX_10GHZ = sp_flux[2]
        EBX_10GHZ = sp_flux[3]
       
        FGL = sp_flux[4]*1000
        EGL = sp_flux[5]*1000

        self.nu_arr = freq_arr
        self.s_nu= flux_arr
        self.e_s_nu = eflux_arr
        self.tmp = freq_arr
        RchiPl = 0
        RchiCpl = 0
        #adding spectral indices.

        
        #observed data points
        # Plotting with errorbars. 
        
        '''
        Making the arrays for the minimization function. This is a piecewise function 
        with spectral indices  specified at negative of the central freq.
        '''
        mymsk = tuple([freq_arr>0])
        ax.errorbar(self.nu_arr[mymsk],self.s_nu[mymsk], yerr=self.e_s_nu[mymsk], 
                    marker ='o', linestyle='none')
        if ALMA[0] >0:
            ax.errorbar(345, ALMA[0], yerr=ALMA[1], marker='d', 
                        color = 'red', linestyle='none')
        #Plotting the bowties.
        if len(alpha_AX)>0 and alpha_AX[1]>-9999:
            if sp_class[ii]['SNR_AX'] > 60:
                rmfit.plot_bowtie(10,alpha_AX, FAX_10GHZ, EAX_10GHZ,ax)
 
        if len(alpha_BX)>0 and alpha_BX[1]>-9999:
            if sp_class[ii]['SNR_BX'] > 60:
                rmfit.plot_bowtie(10,alpha_BX, FBX_10GHZ, EBX_10GHZ,ax)
 
        if len(alpha_GL)>0:
            rmfit.plot_bowtie(0.2005,alpha_GL, FGL,EGL, ax)
            
        # Adjusting the labels, axes scales and range.  
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(0.01*min(flux_arr), 100*max(flux_arr))
        ax.set_xlim(0.001, 1000)
        ax.set_xlabel(r'log ($\nu / GHz)$', fontdict=dict(fontsize=10))
        ax.set_ylabel(r'log ($S_{\nu}/mJy)$', fontdict=dict(fontsize=10))
        ax.tick_params(axis = 'both', which = 'minor', direction='in', length=4, top=True, right=True)
        ax.tick_params(axis = 'both', which = 'major', direction='in', length=9, top=True, right=True)
        ax.text(.6,.8,self.name+'\n'+str(ii)+'\n z ='+str(atgcat['redshift'][ii]),
                horizontalalignment='center',transform=ax.transAxes) #Title

        
        #Initial fitting and plot
        #arrays for fitting:
         
            
        #Guess parameters::: Starting with a guess value. 
        '''
        alpha=-0.7, is the standard value for most radio sources. 
        nu_t, s0: lowest freq available in the observations and its flux. 
        '''
        '''
        s0 = np.max(flux_arr)
        nu_t = np.min(freq_arr)
        alpha = -0.7
        indg = fitTab.loc[self.name].index
        row_g = list(fitTab[indg])
        guess_par = [s0, alpha,nu_t]
        guess_pl = [s0, alpha]
        guess_cpl = [s0,alpha,0.1]
        guess_effa = guess_par; guess_iffa=guess_par; 
        guess_ssa=guess_par
        '''
        
        models = ['PL', 'CPL', 'EFFA', 'IFFA', 'SSA']
        
        model_lists = {'PL': rmfit.power_law, 'CPL': rmfit.curved_power,
                       'EFFA': rmfit.EFFA_func, 'IFFA': rmfit.IFFA_func,
                       'SSA': rmfit.SSA_func, 'GCV': rmfit.gen_curve_func}
        s0 = np.max(flux_arr)
        nu_t = np.min(freq_arr[freq_arr > 0])
        alpha = -0.7
        guess_cpl = [s0, alpha, nu_t]
        guess_pl = [s0, alpha]
        rowFit = fitRes.loc[sname]
        fitFlag = True
        if rowFit[2]>0:
             guess_pl = rowFit[2:4]
             fitFlag = False
        if rowFit[4]>0:
             guess_cpl = rowFit[4:7]
             fitFlag = False

        # Default Mode will plot only two options: CPL and PL
        try:               
            self.fit_pl, self.perr_pl, RchiPl = rmfit.chi_sq_fit(self.nu_arr, self.s_nu,
                                           self.e_s_nu, guess_pl, 'PL')
        except:
            self.fit_pl = [-999, -999]
            
        try:               
            self.fit_cpl, self.perr_cpl, RchiCpl = rmfit.chi_sq_fit(self.nu_arr, self.s_nu,
                                           self.e_s_nu, guess_cpl, 'CPL')
        except:
            self.fit_cpl = [-999, -999, -999]
    
                   
    
        '''
        Calculating reduced chi square values.
        '''
        print(self.nu_arr, self.s_nu,self.e_s_nu, self.fit_pl, RchiPl)
        chi_pl = rmfit.redchisq_calc(model_lists['PL'], self.nu_arr, self.s_nu, 
                                  self.e_s_nu, self.fit_pl)
        chi_cpl = rmfit.redchisq_calc(model_lists['CPL'], self.nu_arr, self.s_nu, 
                                   self.e_s_nu, self.fit_cpl)
 
    
        
        #Sum of the reduced chi square.
        mystr =str(u"\u03C3").upper()+r' $\chi^2_{Red}$'+':\n'
        mystr += 'PL= {0:.2f}\n'.format(RchiPl)
        mystr += 'CPL= {0:.2f}\n'.format(RchiCpl)
        mystr += '**dof for 3par  = {0:d}'.format(len(self.nu_arr) - 4)
        

        
        # Evaluating the model at observed frquencies and plotting it
        self.s1 = rmfit.power_law(self.nu_mod, *self.fit_pl)
        self.s2 = rmfit.curved_power(self.nu_mod, *self.fit_cpl)

        self.l1, = ax.plot(self.nu_mod, self.s1, 'orange', linestyle = '-', label='PL')
        self.l2, = ax.plot(self.nu_mod, self.s2, 'green', linestyle = '-.', label='Curved PL')

        SpShape = sp_class[ii]['Spectral Shape']
        ax.legend(loc = 'upper right')
        self.chistr = axchi.text(0.1,0.05, mystr, wrap=True , fontsize = 9)

        '''Linking all the buttons and sliders with appropriate functions. This
        part makes it all interactive. '''
        # Calling the update function if the slider values are changed. 
        salpha.on_changed(self.update)
        ss0.on_changed(self.update)
        snu.on_changed(self.update)
        
        button.on_clicked(self.reset) # Reset button
        radio.on_clicked(self.modelfunc) # Model selection button
        ref_button.on_clicked(self.refit) #Button forfitting again. 
        bsave.on_clicked(self.save_values) # save button
        brepl.on_clicked(self.replot) #Replot button
    
        #Fit prameters
        pnames_3 = [r'$S_0$', r'$\alpha$',r'$\nu_p$' ]
        mystring = 'Source Info \n'
        mystring += 'Morphology: ' + str(sp_class[ii]['Morphology'] )+ '\n'
        mystring += 'Spectral Shape: ' + str(sp_class[ii]['Spectral Shape'] )+ '\n'
        mystring += 'SNR AX: ' + str(sp_class[ii]['SNR_AX'] )+'\n'
        mystring += 'SNR BX: ' + str(sp_class[ii]['SNR_BX']) +'\n'
        mystring += 'Fitted Parameters \n'        
        mystring += rmfit.par_string(self.fit_pl, self.perr_pl,pnames_3[:2], 'PL' )
        mystring += rmfit.par_string(self.fit_cpl, self.perr_cpl,pnames_3, 'CPL' )
       
        if fig_flag:        
            '''Plotting reduced chi square for the individual data point in the upper panel.'''
            chi_ax.plot(self.tmp, chi_pl, marker='d',color='orange',linestyle='none')
            chi_ax.plot(self.tmp, chi_cpl, marker='*',color='green',linestyle='none')
            chi_ax.tick_params(axis = 'both', which = 'minor', direction='in', 
                               length=4, top=True, right=True)
            chi_ax.tick_params(axis = 'both', which = 'major', direction='in', 
                               length=9, top=True, right=True)
            chi_ax.set_xscale('log')
            chi_ax.set_yscale('log')      
            chi_ax.set_xlim(0.01*min(freq_arr), 100*max(freq_arr))
            chi_ax.set_ylabel(r'Red $\chi^2$')
            
            ###################################################################
            
            # This part plots two point alpha in the alpha axis. 
            freq_alph = []
            alpha_arr  = []
            ealpha_arr = []
            
            alpha_tab = Table(data=(freq_arr, flux_arr, eflux_arr), names=('Freq', 'Flux', 'err_flux'))
            alpha_tab.sort('Freq')
            for ial in range(0,len(alpha_tab)-1):
                nu_al = [alpha_tab['Freq'][ial], alpha_tab['Freq'][ial+1]]
                if alpha_tab['Freq'][ial] != alpha_tab['Freq'][ial+1]:                
                    s_al = [alpha_tab['Flux'][ial], alpha_tab['Flux'][ial+1]]
                    err_al = [alpha_tab['err_flux'][ial], alpha_tab['err_flux'][ial+1]]
                    mid_al = rmfit.alpha_calc(nu_al, s_al, err_al)
                    freq_alph.append(np.sqrt(nu_al[0]*nu_al[1]))
                    alpha_arr.append(mid_al[0])
                    ealpha_arr.append(mid_al[1])
                    
            alpha_ax.errorbar(freq_alph,alpha_arr, yerr=ealpha_arr, linestyle='none', 
                              marker='d')
            
            if len(alpha_AX)>0:
                if sp_class[ii]['SNR_AX'] > 60:
                    alpha_ax.errorbar(10,alpha_AX[1], yerr=alpha_AX[2], linestyle='none', marker='*', 
                                  color='red', label='JVLA AX')
            if len(alpha_BX)>0:
                if sp_class[ii]['SNR_BX'] > 60:
                    alpha_ax.errorbar(10,alpha_BX[1], yerr=alpha_BX[2], linestyle='none', marker='o',
                                  color='green', label= 'JVLA BX')     
            if len(alpha_GL)>0:
                alpha_ax.errorbar(0.2005,alpha_GL[1], yerr=alpha_GL[2], linestyle='none', marker='d', 
                                  color='orange', label='GLEAM')
                
            if plot_flag == 'all':    
                alpha_ax.plot(self.fit_effa[2], 0, '^', color='blue', alpha=0.5)
                alpha_ax.plot(self.fit_iffa[2], 0, '^', color='red',alpha=0.5)
                alpha_ax.plot(self.fit_ssa[2], 0, '^', color='green',alpha=0.5)
            alpha_ax.plot(self.fit_cpl[2], 0, '^', color='blue', alpha=0.5)        
            alpha_ax.legend()

        #######################################################################
        #######################################################################
        '''######## A bit to plot nu_p from CPL:'''
        ########
        ''' x = 1/2(-alpha/q +/- sqrt(lnbeta/q))  where x = ln(nu_+/-'''
        '''######## A bit to plot nu_p from CPL:'''
        
        cpl_q = self.fit_cpl[2]
        cpl_a = self.fit_cpl[1]
        cpl_s = self.fit_cpl[0]
        cpl_q_err = self.perr_cpl[2]
        cpl_a_err = self.perr_cpl[1]
        cpl_s_err = self.perr_cpl[0]
        nu_to = -999
        nu_to_err = -999
        beta = -999
        ########
        ''' x = 1/2(-alpha/q +/- sqrt(lnbeta/q))  where x = ln(nu_+/-'''

        if SpShape =='CV' or SpShape=='CV:':
            nu_p = np.exp(-cpl_a/(2*cpl_q))
            fvar_q = (cpl_q_err/cpl_q)**2
            fvar_a = (cpl_a_err/cpl_a)**2
            fvar_s = (cpl_s_err/cpl_s)**2
            f_dec = np.sqrt(fvar_q + fvar_a + fvar_s)
            beta = 1 - f_dec
            fac = np.exp(np.sqrt(np.log(beta)/cpl_q))
            vlim_1 = nu_p*fac 
            vlim_2 = nu_p/fac
            #vlim_1 = np.exp(0.5*((-cpl_a/cpl_q) + np.sqrt(np.log(beta)/cpl_q)))
            #vlim_2 = np.exp(0.5*((-cpl_a/cpl_q) - np.sqrt(np.log(beta)/cpl_q)))
            l_lim = nu_p - min(vlim_1, vlim_2)
            u_lim = max(vlim_1, vlim_2) - nu_p

            s_nu_p = rmfit.curved_power(nu_p, *self.fit_cpl)
            s_nu_p_err = s_nu_p*f_dec
            ax.errorbar(nu_p, s_nu_p, xerr=[[l_lim], [u_lim]],yerr=s_nu_p_err, linestyle='none', marker = 'd', 
                     color='black',  fillstyle='none')
            mystring += r'$\beta$ = ' +str(np.round(beta,2)) + '\n'
            mystring += r'$\nu_p$ = ' + str(np.round(nu_p,2)) + r'$\pm$' +str(np.round(l_lim,2)) +'\n'
            nu_to=nu_p
            nu_to_err = l_lim
               
        
        elif SpShape=='PK' or SpShape=='PK:':
            nu_p = np.exp(-cpl_a/(2*cpl_q))
            fvar_q = (cpl_q_err/cpl_q)**2
            fvar_a = (cpl_a_err/cpl_a)**2
            fvar_s = (cpl_s_err/cpl_s)**2
            f_dec = 0.5*np.sqrt(fvar_q + fvar_a + fvar_s)
            beta = 1 - f_dec
            fac = np.exp(np.sqrt(np.log(beta)/cpl_q))
            vlim_1 = nu_p*fac
            vlim_2 = nu_p/fac

            #vlim_1 = np.exp(0.5*((-cpl_a/cpl_q) + np.sqrt(np.log(beta)/cpl_q)))
            #vlim_2 = np.exp(0.5*((-cpl_a/cpl_q) - np.sqrt(np.log(beta)/cpl_q)))
            l_lim = nu_p - min(vlim_1, vlim_2)
            u_lim = max(vlim_1, vlim_2) - nu_p
            s_nu_p = rmfit.curved_power(nu_p, *self.fit_cpl)
            s_nu_p_err = s_nu_p*f_dec
            ax.errorbar(nu_p, s_nu_p, xerr=[[l_lim], [u_lim]],yerr=s_nu_p_err, linestyle='none', marker = 'd', 
                     color='black',  fillstyle='none')
            mystring += r'$\beta$ = ' +str(np.round(beta,2)) + '\n'
            mystring += r'$\nu_p$ = ' + str(np.round(nu_p,2)) + r'$\pm$' +str(np.round(l_lim,2)) +'\n'
            nu_to=nu_p
            nu_to_err = l_lim
            
        elif SpShape=='PL' or SpShape=='PL:':
            onu = np.array([freq_arr, flux_arr])
            snu_arr = np.sort(onu)
            nu_p = snu_arr[0,0]/10.0
            s_nu_p = rmfit.power_law(nu_p, *self.fit_pl)
            ax.errorbar(nu_p, s_nu_p, xerr=nu_p*0.5, linestyle='none', marker = 'd', 
                     color='black', xuplims = True, fillstyle='none')
            nu_to=nu_p

        
        #######################################################################
        #######################################################################
        if fig_flag: 
            '''
            FITSfigure
        
            '''
            image_width = 8.0/3600
    
            zscale = ZScaleInterval()
            source_name = rmfit.func_iau_names(self.name)
            fits_ax = fitsdir_AX+source_name+'F.AX.fits'
            fits_bx = fitsdir_BX+source_name+'F.BX.fits'
            
            sflag_ax = len(ax_list[ax_list.isin([source_name])])
            sflag_bx = len(bx_list[bx_list.isin([source_name])])
            if sflag_ax == 0:
                box = plt.axes([0.77,.35,0.22,0.22], frameon = True )
                box.tick_params(axis='both', bottom='off', left='off', labelbottom='off', 
                                labelleft='off')
                box.text( 0.2, 0.5,'The target was not observed.', wrap = True)
                box.set_title(source_name+'.AX')
    
            
            elif os.path.isfile(fits_ax):
                if len(jvla_AX)>1:
                    jvla_AX = jvla_AX[0]
                f1 = apl.FITSFigure(fits_ax, figure = fig, subplot= [0.77,.35,0.22,0.22])
                hdu = fits.open(fits_ax)
                data = hdu[0].data
                headers = hdu[0].header
                hdu.close()
                xpm = headers['NAXIS1']/2
                ypm = headers['NAXIS2']/2 
                rac, decc = f1.pixel2world(xpm,ypm)
                ax_cord = SkyCoord(rac, decc, frame = 'fk5',unit=(u.deg, u.deg))            
                ra = ax_cord.ra.value
                dec = ax_cord.dec.value
                dra, ddec = rmfit.tick_offsets(f1, ax_cord)
                zmin, zmax = zscale.get_limits(data)
                f1.show_colorscale(cmap=plt.viridis(), vmin = zmin, vmax = zmax, stretch='power', exponent = 2)
                f1.recenter(ra,dec, width = image_width, height = image_width )
                f1._ax1.xaxis.set_ticklabels(dra)
                f1._ax1.yaxis.set_ticklabels(ddec)
                f1.axis_labels.set_xtext(r'$\Delta$'+' RA asec' )
                f1.axis_labels.set_ytext(r'$\Delta$'+' DEC asec' )
                f1.set_title(source_name+'.AX')
                f1.add_beam()
                f1.beam.set_corner('bottom left')
                f1.beam.set_color('red')
                f1.beam.set_frame(True)
                f1.show_markers(atgcat['wise_ra'][ii], atgcat['wise_dec'][ii])
                f1.ticks.set_minor_frequency(5)
                f1.ticks.set_length(5)
                f1.ticks.set_color('white')
                f1._ax1.tick_params(axis='both',which='both',  direction='in')
                f1._ax2.tick_params(axis='both',which='both', direction='in')
                f1.ticks.set_linewidth(2)
                
            else:
                box = plt.axes([0.77,.35,0.22,0.22], frameon = True )
                box.tick_params(axis='both', bottom='off', left='off', labelbottom='off', 
                                labelleft='off')
                box.text( 0.2, 0.5,'The target image has some issues.', wrap = True)
                box.set_title(source_name+'.AX')
    
    
            if sflag_bx == 0:
                box = plt.axes([0.77,0.68,0.22,0.22], frameon = True )
                box.tick_params(axis='both', bottom='off', left='off', labelbottom='off', 
                                labelleft='off')
                box.text( 0.2, 0.5,'The target was not observed', wrap = True)
                box.set_title(source_name+'.BX')
    
            elif os.path.isfile(fits_bx):
                if len(jvla_BX)>1:
                    jvla_BX = jvla_BX[0]
                f2 = apl.FITSFigure(fits_bx, figure = fig, subplot= [0.77,0.68,0.22,0.22])
                hdu = fits.open(fits_bx)
                data = hdu[0].data
                headers = hdu[0].header
                hdu.close()
                xpm = headers['NAXIS1']/2
                ypm = headers['NAXIS2']/2 
                rac, decc = f2.pixel2world(xpm,ypm)
                bx_cord = SkyCoord(rac, decc, frame = 'fk5',unit=(u.deg, u.deg))            
                ra = bx_cord.ra.value
                dec = bx_cord.dec.value
                dra_b, ddec_b = rmfit.tick_offsets(f2, bx_cord)
                zmin, zmax = zscale.get_limits(data)
                f2.show_colorscale(cmap=plt.viridis(), vmin = zmin, vmax = zmax, 
                                   stretch='power', exponent = 2)
                f2.recenter(ra,dec, width = image_width, height = image_width )
                f2._ax1.xaxis.set_ticklabels(dra_b)
                f2._ax1.yaxis.set_ticklabels(ddec_b)
                f2.axis_labels.set_xtext(r'$\Delta$'+' RA asec' )
                f2.axis_labels.set_ytext(r'$\Delta$'+' DEC asec' )
                f2.set_title(source_name+'.BX')
                f2.add_beam()
                f2.beam.set_corner('bottom left')
                f2.beam.set_color('red')
                f2.beam.set_frame(True)
                f2.show_markers(atgcat['wise_ra'][ii], atgcat['wise_dec'][ii])
                f2.ticks.set_minor_frequency(5)
                f2.ticks.set_length(5)
                f2.ticks.set_color('white')
                f2._ax1.tick_params(axis='both',which='both',  direction='in')
                f2._ax2.tick_params(axis='both',which='both', direction='in')
                f2.ticks.set_linewidth(2)
                
            else:
                box = plt.axes([0.77,0.68,0.22,0.22], frameon = True )
                box.tick_params(axis='both', bottom='off', left='off', labelbottom='off', 
                                labelleft='off')
                box.text( 0.2, 0.5,'The target image has some issues.', wrap = True)
                box.set_title(source_name+'.BX')
                
                
            
        #######################################################################
        #######################################################################
        ''' Extrapolation at ALMA:'''
        nu_ALMA = 345 # In GHZ
        S_ALMA = -999
        S_ALMA_err = -999
        f_syn = -999
        f_syn_err = -999
        alpha_345 = -999
        alpha_345_err = -999
        if SpShape =='CV' or SpShape=='CV:' or SpShape=='PK' or SpShape=='PK:':
            
            #identify unique two freq and flux points:
            nu_mm = [10]
            if FAX_10GHZ>0:
                snu_mm = [FAX_10GHZ]
                esnu_mm = [EAX_10GHZ]
            else:
                snu_mm = [FBX_10GHZ]
                esnu_mm = [EBX_10GHZ]
                
            low_nu =np.sort(np.unique(self.nu_arr))[-2]
            low_ind = np.argwhere(self.nu_arr ==low_nu)
            nu_mm.append(low_nu)
            snu_mm.append(self.s_nu[low_ind])
            esnu_mm.append(self.e_s_nu[low_ind])
            
        
            fit_alma_pl,cov_a = scipy.optimize.curve_fit(rmfit.power_law, 
                        np.array(nu_mm),np.array(snu_mm), self.fit_pl,
                                          np.array(esnu_mm))
           
            perr_alma_pl = np.sqrt(np.diag(cov_a))
            alpha_345 = fit_alma_pl[1]
            alpha_345_err = perr_alma_pl[1]
            
            
            S_ALMA = rmfit.power_law(nu_ALMA, *fit_alma_pl)
            S_ALMA_err = S_ALMA*((perr_alma_pl[0]/fit_alma_pl[0]) +
                             (np.log(nu_ALMA)*perr_alma_pl[1]))
            ax.errorbar(nu_ALMA, S_ALMA, yerr=S_ALMA_err, linestyle='none', marker = '^', 
                 color='black', xuplims = True, fillstyle='none')
            if ALMA[0]>0:
                f_syn = S_ALMA/ALMA[0]
                f_syn_err = f_syn*((S_ALMA_err/S_ALMA) + ALMA[1]/ALMA[0])
                mystring += '\n'+r'f$_{syn}^{345GHz}$ = '
                mystring += str(np.round(f_syn,2)) + r'$\pm$' + str(np.round(f_syn_err, 2))+'\n'
            
        elif SpShape =='PL' or SpShape =='PL:':
            print('Power_Law')
            S_ALMA = rmfit.power_law(nu_ALMA, *self.fit_pl)
            S_ALMA_err = S_ALMA*((self.perr_pl[0]/self.fit_pl[0]) +
                             (np.log(nu_ALMA)*self.perr_pl[1]))
            
            if fig_flag:            
                ax.errorbar(nu_ALMA, S_ALMA, yerr=S_ALMA_err, linestyle='none', marker = '^', 
                     color='black', xuplims = True, fillstyle='none')
            if ALMA[0]>0:
                f_syn = S_ALMA/ALMA[0]
                f_syn_err = f_syn*((S_ALMA_err/S_ALMA) + ALMA[1]/ALMA[0])
                mystring += '\n'+r'f$_{syn}^{345GHz}$ = '
                mystring += str(np.round(f_syn,2)) + r'$\pm$' + str(np.round(f_syn_err, 2))+'\n'

            
        if fig_flag: 
            infoax.text(0.03,0.00,mystring, wrap=True, fontsize=9)
            plt.savefig('./Plots_SNRF/'+self.name+'.png')

            
        #######################################################################
        #######################################################################


    def update(self,val):
        '''If value in the slider is changed manually or through a function,
        replot the given/selected model'''
        alpha = salpha.val
        s0 = 10**ss0.val
        nu_t = 10**snu.val
        if radio.value_selected=='EFFA':
            s_up = rmfit.EFFA_func(self.nu_mod, s0, alpha, nu_t)
            self.l1.set_ydata(s_up)
        if radio.value_selected=='IFFA':
            s_up = rmfit.IFFA_func(self.nu_mod, s0, alpha, nu_t)
            self.l2.set_ydata(s_up)
        if radio.value_selected=='SSA':
            s_up = rmfit.SSA_func(self.nu_mod, s0, alpha, nu_t)
            self.l3.set_ydata(s_up)             
        if radio.value_selected=='PL':
            s_up = rmfit.power_law(self.nu_mod, s0, alpha)
            self.l4.set_ydata(s_up)    
        if radio.value_selected=='CPL':
            s_up = rmfit.curved_power(self.nu_mod, s0, alpha, snu.val)
            self.l5.set_ydata(s_up) 
                        
        fig.canvas.draw_idle()
    
    def reset(self,event):
        '''Reset to the starting value'''
        salpha.reset()
        ss0.reset()
        snu.reset()
        if radio.value_selected=='PL':
            self.l4.set_ydata(self.s4)
        if radio.value_selected=='CPL':
            self.l4.set_ydata(self.s5)            
        if radio.value_selected=='EFFA':
            self.l1.set_ydata(self.s1)
        if radio.value_selected=='IFFA':
            self.l2.set_ydata(self.s2)
        if radio.value_selected=='SSA':
            self.l3.set_ydata(self.s3)           
        fig.canvas.draw_idle()

     
    def modelfunc(self,label):
        '''This allows to play with only one model at at time.'''
        if label=='EFFA':
            s_up = rmfit.EFFA_func(self.nu_mod, self.fit_effa)
            self.l1.set_ydata(s_up)
            self.l1.set_color('blue')
            self.l2.set_visible(False)
            self.l3.set_visible(False)
            self.l4.set_visible(False)
            self.l5.set_visible(False)

        if label=='IFFA':
            s_up = rmfit.IFFA_func(self.nu_mod, self.fit_iffa)
            self.l2.set_ydata(s_up)
            self.l2.set_color('red')
            self.l1.set_visible(False)
            self.l3.set_visible(False)
            self.l4.set_visible(False)
            self.l5.set_visible(False)

        if label=='SSA':
            s_up = rmfit.SSA_func(self.nu_mod, self.fit_ssa)
            self.l3.set_ydata(s_up)
            self.l3.set_color('green')
            self.l2.set_visible(False)
            self.l1.set_visible(False)
            self.l4.set_visible(False)
            self.l5.set_visible(False)
        if label=='PL':           
            s_up = rmfit.power_law(self.nu_mod, self.fit_pl)
            self.l4.set_ydata(s_up)
            self.l4.set_color('orange')
            self.l5.set_visible(False)
            self.l2.set_visible(False)
            self.l3.set_visible(False)
            self.l1.set_visible(False)
        if label=='CPL':           
            s_up = rmfit.curved_power(self.nu_mod, self.fit_cpl)
            self.l5.set_ydata(s_up)
            self.l5.set_color('green')
            self.l4.set_visible(False)
            self.l2.set_visible(False)
            self.l3.set_visible(False)
            self.l1.set_visible(False)
        

        fig.canvas.draw_idle()
        
        
        
        
    def refit(self,event):
        '''If new guess values are looking better, refit function will run the 
        curve_fit again, and print out the new best estimates.  Need to press
        c key to start the fit. Till then one can play with the replot.'''
        label=radio.value_selected
        print (label+'\n')
        print('Refitting with new guess values\n')
        
        print("Press c to fit with the new guess parameters.")
        print("Press q to stop")
        presser = fig.canvas.mpl_connect('key_press_event',
                                             self.onkeypress)       
        while True:
            # this function returns True if a key was pressed,
            # False if a mouse button was clicked, and None
            # if neither happened within timeout
            keypressed = fig.waitforbuttonpress()
            '''If only c key is pressed do the following'''
            if keypressed and self.keypress =='c':
                nu_peak= 10**snu.val
                s0 = 10**ss0.val
                alpha = salpha.val
                label=radio.value_selected
                guess_par=[s0, alpha, nu_peak]
                guess_pl = [s0, alpha]
                guess_cpl= [s0, alpha, snu.val]
                if label=='EFFA':
                    self.fit_effa,cov1 = scipy.optimize.curve_fit(rmfit.model_effa, self.nu_arr, self.s_nu, guess_par, self.e_s_nu)   
                    self.perr_effa = np.sqrt(np.diag(cov1))
                    ss0.set_val(np.log10(self.fit_effa[0]))
                    snu.set_val(np.log10(self.fit_effa[2]))
                    salpha.set_val(self.fit_effa[1])
                    print ('Func: '+label+'\n'+'Guess Values: ')
                    print(guess_par)
                if label=='IFFA':
                    self.fit_iffa,cov2 = scipy.optimize.curve_fit(rmfit.model_iffa, self.nu_arr, self.s_nu, guess_par, self.e_s_nu)        
                    self.perr_iffa = np.sqrt(np.diag(cov2))
                    ss0.set_val(np.log10(self.fit_iffa[0]))
                    snu.set_val(np.log10(self.fit_iffa[2]))
                    salpha.set_val(self.fit_iffa[1])
                    print ('Func: '+label+'\n'+'Guess Values: ')
                    print(guess_par)
                if label=='SSA':
                    self.fit_ssa, cov3 = scipy.optimize.curve_fit(rmfit.model_ssa,  self.nu_arr, self.s_nu, guess_par, self.e_s_nu)        
                    self.perr_ssa = np.sqrt(np.diag(cov3))
                    ss0.set_val(np.log10(self.fit_ssa[0])) 
                    snu.set_val(np.log10(self.fit_ssa[2]))
                    salpha.set_val(self.fit_ssa[1])
                    print ('Func: '+label+'\n'+'Guess Values: ')
                    print(guess_par)
                if label=='PL':
                    self.fit_pl, cov4 = scipy.optimize.curve_fit(rmfit.model_pl,  self.nu_arr, self.s_nu, guess_pl, self.e_s_nu)        
                    self.perr_pl = np.sqrt(np.diag(cov4))
                    ss0.set_val(np.log10(self.fit_pl[0])) 
                    #snu.set_val(np.log10(self.fit_pl[2]))
                    salpha.set_val(self.fit_pl[1])
                    print ('Func: '+label+'\n'+'Guess Values: ')
                    print([s0,alpha])
                    
                if label=='CPL':
                    self.fit_cpl, cov5 = scipy.optimize.curve_fit(rmfit.model_cpl,  self.nu_arr, self.s_nu, guess_cpl, self.e_s_nu)        
                    self.perr_cpl = np.sqrt(np.diag(cov5))
                    ss0.set_val(np.log10(self.fit_cpl[0])) 
                    snu.set_val(self.fit_cpl[2])
                    salpha.set_val(self.fit_cpl[1])
                    print ('Func: '+label+'\n'+'Guess Values: ')
                    print([s0,alpha,snu.val])

                   
                
                if keypressed  and self.keypress =='q' :
                    break
                '''
                Printing new values.
                '''
                pnames_3 = ['$S_0$', '$\alpha$','$\nu_p$' ]
                pnames_cpl = [r'$S_0$', r'$\alpha$','q' ]
                mystring = 'Fitted Parameters \n'
                mystring += rmfit.par_string(self.fit_effa, self.perr_effa, pnames_3, 'EFFA')
                mystring += rmfit.par_string(self.fit_iffa, self.perr_iffa, pnames_3, 'IFFA')
                mystring += rmfit.par_string(self.fit_ssa, self.perr_ssa, pnames_3, 'SSA')
                mystring += rmfit.par_string(self.fit_pl, self.perr_pl,pnames_3[:2], 'PL' )
                mystring += rmfit.par_string(self.fit_cpl, self.perr_cpl,pnames_cpl, 'CPL' )

                infoax.clear()
                infoax.text(0.03,0.00,mystring, wrap=True, fontsize=9 ) 
                chi_ax.clear()
                #print (self.nu_arr, self.s_nu, self.e_s_nu)
                '''
                Calculating reduced chi square values.
                '''
                chi_effa = rmfit.redchisq_calc(rmfit.model_effa,self.nu_arr, self.s_nu, self.e_s_nu, self.fit_effa )
                chi_iffa = rmfit.redchisq_calc(rmfit.model_iffa,self.nu_arr, self.s_nu, self.e_s_nu, self.fit_iffa )
                chi_ssa = rmfit.redchisq_calc(rmfit.model_ssa,self.nu_arr, self.s_nu, self.e_s_nu, self.fit_ssa )
                chi_pl = rmfit.redchisq_calc(rmfit.model_pl, self.nu_arr, self.s_nu, self.e_s_nu, self.fit_pl)
                chi_cpl = rmfit.redchisq_calc(rmfit.model_cpl, self.nu_arr, self.s_nu, self.e_s_nu, self.fit_cpl)
         
                #Sum of the reduced chi square.
                mystr =str(u"\u03C3").upper()+r' $\chi^2_{Red}$'+':\n'
                mystr += 'EFFA= '+str(np.round(np.sum(chi_effa),2))+'\n'
                mystr += 'IFFA= '+str(np.round(np.sum(chi_iffa),2))+'\n'
                mystr += 'SSA= '+str(np.round(np.sum(chi_ssa),2))+'\n'
                mystr += 'PL= '+str(np.round(np.sum(chi_pl),2))+'\n'
                mystr += 'CPL= '+str(np.round(np.sum(chi_cpl),2))+'\n'
                mystr += '**dof for 3par  = '+str(len(self.nu_arr) - 4)

                axchi.clear()
                self.chistr = axchi.text(0.1,0.05, mystr, wrap=True , fontsize = 9)
                
                chi_ax.plot(self.tmp, chi_effa,marker='o',color='blue',linestyle='none')
                chi_ax.plot(self.tmp, chi_iffa,marker = 's',color='red',linestyle='none')
                chi_ax.plot(self.tmp, chi_ssa, marker='^',color='green',linestyle='none')
                chi_ax.plot(self.tmp, chi_pl, marker='d',color='orange',linestyle='none')
                chi_ax.plot(self.tmp, chi_cpl, marker='*',color='green',linestyle='none')

                chi_ax.tick_params(axis = 'both', which = 'minor', direction='in', length=4, top=True, right=True)
                chi_ax.tick_params(axis = 'both', which = 'major', direction='in', length=9, top=True, right=True)
                chi_ax.set_xscale('log')
                chi_ax.set_yscale('log')        
                chi_ax.set_xlim(0.01*min(np.abs(self.nu_arr)), 100*max(np.abs(self.nu_arr)))
                chi_ax.set_ylabel(r'Red $\chi^2$')

                print (mystring)
                fig.canvas.draw_idle()

        # kill the event watchers
        #fig.canvas.mpl_disconnect(clicker)
        fig.canvas.mpl_disconnect(presser)
        
        
    def replot(self,event):
        '''Plot the models using selected parameters.'''
        label=radio.value_selected
        print (label+'\n')
        print('Replot with new guess values\n')
        
        print("Press r to replot the model.")
        print("Press q to stop")
        # set up the key-press and mouse-click event watcher
        #clicker = fig.canvas.mpl_connect('button_press_event',
        #                                      self.onclick)       
        presser = fig.canvas.mpl_connect('key_press_event',
                                             self.onkeypress)       
        while True:
            # this function returns True if a key was pressed,
            # False if a mouse button was clicked, and None
            # if neither happened within timeout
            keypressed = fig.waitforbuttonpress()
            '''Press r to select the turnover and use slider for alpha'''
            if keypressed and self.keypress =='r':
                nu_peak= self.clickx_data
                S_nu_peak = self.clicky_data
                label=radio.value_selected
                if label=='EFFA':
                    alpha = self.fit_effa[1]
                    s0= S_nu_peak/(nu_peak**alpha*np.exp(-1.0))
                    ss0.set_val(np.log10(s0))
                    snu.set_val(np.log10(nu_peak))
                    
                if label=='IFFA':
                    alpha = self.fit_iffa[1]                   
                    s0 = S_nu_peak/(nu_peak**alpha*(1-np.exp(-(1.0)**(-2.1))))
                    ss0.set_val(np.log10(s0))
                    snu.set_val(np.log10(nu_peak))


                if label=='SSA':
                    tau=1
                    alpha = self.fit_ssa[1]
                    s0= S_nu_peak/(1-np.exp(-tau))
                    ss0.set_val(np.log10(s0))
                    snu.set_val(np.log10(nu_peak))
                     
                if label=='PL':
                    s0= S_nu_peak
                    nu_0 = nu_peak
                    ss0.set_val(np.log10(s0))
                    snu.set_val(np.log10(nu_0))
                    
                if label=='CPL':
                    s0= S_nu_peak
                    nu_0 = nu_peak
                    ss0.set_val(np.log10(s0))
                    snu.set_val(nu_0)

                   
                if keypressed  and self.keypress =='q' :
                    break
                
                fig.canvas.draw_idle()
                

        # kill the event watchers
        #fig.canvas.mpl_disconnect(clicker)
        fig.canvas.mpl_disconnect(presser)




        #Fit prameters
        
    def save_values(self, event):
        label=radio.value_selected
        ind = fitTab.loc[self.name].index

        if label=='EFFA':
            fitTab[ind]['EFFA_s0']=self.fit_effa[0]
            fitTab[ind]['EFFA_alpha']=self.fit_effa[1]
            fitTab[ind]['EFFA_nup']=self.fit_effa[2]
        if label=='IFFA':
            fitTab[ind]['IFFA_s0']=self.fit_iffa[0]
            fitTab[ind]['IFFA_alpha']=self.fit_iffa[1]
            fitTab[ind]['IFFA_nup']=self.fit_iffa[2]
        if label=='SSA':
            fitTab[ind]['SSA_s0']=self.fit_ssa[0]
            fitTab[ind]['SSA_alpha']=self.fit_ssa[1]
            fitTab[ind]['SSA_nup']=self.fit_ssa[2]
        ascii.write(fitTab, output='GuessPar_Radiofits.csv', format='csv', overwrite=True)

        
        
    '''example functions to record keypress and click events'''    
    def onclick(self,event):
        """
        Handle mouse click event
        """
        if not event.inaxes:
            # skip if event happens outside plot
            return
        print("Mouse button {0} pressed at ({1:.2f},{2:.2f})".format(event.button,event.xdata,event.ydata))
        #self.keypress.append(None) # no key press
        #self.clickbutton.append(event.button)
        #self.clickx_data.append(event.xdata)
        #self.clicky_data.append(event.ydata)

    def onkeypress(self,event):
        """
        Handle key press event
        """
        if not event.inaxes:
            # skip if event happens outside plot
            return
        print("Key {0} pressed at ({1:.2f},{2:.2f})".format(event.key,event.xdata,event.ydata))
        self.keypress = event.key
        #self.clickbutton.append(None) # no mouse click
        self.clickx_data = event.xdata
        self.clicky_data = event.ydata
        
        
        
 
    def to_selection(self,event):
        """
        Turnover selection. 
        """
        print("Press r in the plot to selecr turnover frequency .")
        print("Press q to stop")
        # set up the key-press and mouse-click event watcher
        clicker = fig.canvas.mpl_connect('button_press_event',
                                              self.onclick)       
        presser = fig.canvas.mpl_connect('key_press_event',
                                             self.onkeypress)       
        while True:
            # this function returns True if a key was pressed,
            # False if a mouse button was clicked, and None
            # if neither happened within timeout
            keypressed = fig.waitforbuttonpress()
            if keypressed and self.keypress =='r':
                nu_peak= self.clickx_data
                S_nu_peak = self.clicky_data
                alpha = salpha.val
                label=radio.value_selected
                if label=='EFFA':
                    s0= S_nu_peak/(nu_peak**alpha*np.exp(-1.0))
                if label=='IFFA':
                    s0 = S_nu_peak/(nu_peak**alpha*(1-np.exp(-(1.0)**(-2.1))))
                if label=='SSA':
                    tau=1
                    s0= S_nu_peak/(1-np.exp(-tau))
                if label=='PL':
                    s0= S_nu_peak
                
                self.refit([s0, alpha, nu_peak])
            # stop if a key is pressed and that key is "q"
            if keypressed  and self.keypress =='q' :
                break
        # kill the event watchers
        fig.canvas.mpl_disconnect(clicker)
        fig.canvas.mpl_disconnect(presser)
    '''If clicked on next or previous button, clear the figure and run the plotting'''
    def next(self, event):
        self.ind += 1
        i = self.ind
        ax.clear()
        infoax.clear()
        chi_ax.clear()
        axchi.clear()
        alpha_ax.clear()
        self.plotting(i)
    def prev(self, event):
        self.ind -= 1
        i = self.ind 
        ax.clear()
        infoax.clear()
        chi_ax.clear()
        axchi.clear()
        alpha_ax.clear()
        self.plotting(i)
        
    def fastforward(self,event):
        self.ind += 10
        i = self.ind
        ax.clear()
        infoax.clear()
        chi_ax.clear()
        axchi.clear()
        alpha_ax.clear()
        self.plotting(i)

    def fastbackward(self,event):
        self.ind -= 10
        i = self.ind
        ax.clear()
        infoax.clear()
        chi_ax.clear()
        axchi.clear()
        alpha_ax.clear()
        self.plotting(i)


i = 0
obj_pl = Plot_class()

obj_pl.plotting(i)
bnext.on_clicked(obj_pl.next)
bprev.on_clicked(obj_pl.prev)
fastfor.on_clicked(obj_pl.fastforward)
fastbck.on_clicked(obj_pl.fastbackward)
