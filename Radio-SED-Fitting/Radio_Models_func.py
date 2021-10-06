#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 10:33:07 2018

@author: pallavipatil
"""

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import polyfit
from scipy.optimize import curve_fit


class RadioModelFit:
    # This file includes definitions of all radio spectral emission models
    ##########################################################################
    # Adding power law and sign-of curvature::
    def __init__(self):
        self.fit_res = []
        self.perr_res = []
        self.rchiq = None
        self.model_name = 'PL'

    def power_law(self, nu, s0, alpha):
        return s0 * (nu) ** alpha

    def curved_power(self, nu, s0, alpha, nu_t):
        q = alpha / (2 * np.log(nu_t))
        return s0 * (nu ** alpha) * np.exp(q * np.log(nu) * np.log(nu))
    
    def CPL(self, nu, s0, alpha, q):
        return s0 * (nu ** alpha) * np.exp(q * np.log(nu) * np.log(nu))
 

    ###########################################################################

    # External Free Free Absorption:
    # A screen of ionised gas; external to the gas system can lead to the
    # absorption.
    # S_nu = s0* nu**-alpha * exp-[(nu/nu_t)**(-2.1)]
    # nu_t - is the turnover frequency where the optical depth is unity.
    # s0 is the normalization factor
    # alpha is the index of the non-thermal power low

    def EFFA_func(self, nu, s0, alpha, nu_t):
        return s0 * nu ** alpha * np.exp(-(nu / nu_t) ** (-2.1))

    ##########################################################################
    # Internal Free Free Absorption:
    # Thermal plasma mixed with relativistic particles can produce spectral TO.
    # Similar to escape probability mechanism. 
    # S_nu = s0* nu**-alpha * [(1-exp(-tau)) / (tau)]
    # Where tau ~ nu**(-2.1) 

    def IFFA_func(self, nu, s0, alpha, nu_t):
        tau = (nu/nu_t)**(-2.1)
        return s0 * nu ** alpha * ((1 - np.exp(-tau))/tau)
                

    ##########################################################################
    # Synchrotron Self Absorption
    # The source TO is becuase the source cannot have a brightness temp greater
    # than the plasma temperature of the nonthermal electrons. 
    # alpha = -(beta-1)/2--- beta = 1-2*alpha 
    # S_nu = s0 *(nu/nu_t)**-(beta-1)/2 * (1-exp**(-tau)/ tau)
    # Where tau = (nu/nu_t)**-(beta+4)/2

    def SSA_func(self, nu, s0, alpha, nu_t):
        beta = 1 - 2 * alpha
        tau = (nu / nu_t) ** (-(beta + 4) / 2)
        return s0 * ((nu / nu_t) ** (-(beta - 1) / 2)) * (
                (1 - np.exp(-tau)) / tau)

    def gen_curve_func(self, nu, s0, athick, athin, nu_t):
        amp = s0 / (1 - (1 / np.e))
        thick = (nu / nu_t) ** athick
        expo = (nu / nu_t) ** (athin - athick)
        thin = 1 - np.exp(-expo)
        return amp * thin * thick

    def spec_index_calc(self, x, y, err):
        X = np.log10(np.array(x))
        Y = np.log10(np.array(y))
        sol = polyfit(X, Y, 1)
        alpha = sol[0]
        sigma_alpha = np.sqrt((err[0] ** 2) * (
                1.0 / (np.log10(x[0] / x[1])) * 1.0 / (
                np.log(10) * y[0])) ** 2 + (err[1] ** 2) * (1.0 / (
            np.log10(x[0] / x[1])) * -1.0 / (np.log(10) * y[1])) ** 2)
        return alpha, sigma_alpha

    def plot_bowtie(self, nu, alpha_x, flux, e_flux, ax):
        """ A bowtie plotting function. This is used to show in-band spectral
        indices
        and their errors on the plot. """
        s0 = flux
        nu_t = abs(alpha_x[0])
        sp1 = alpha_x[1]
        sp2 = sp1 + alpha_x[2]
        sp3 = sp1 - alpha_x[2]
        # nu_cen = np.linspace(nu/2,2*nu,num=5)
        nu_cen = np.linspace(nu - 2, nu + 2, num=5)
        # f1 = s0*np.power((nu_cen/nu_t), sp1)
        f2 = s0 * np.power((nu_cen / nu_t), sp2)
        f3 = s0 * np.power((nu_cen / nu_t), sp3)
        # ax.plot(nu_cen,f1,'-', linewidth = 1.5)
        ax.plot(nu_cen, f2, '-', color='k', )
        ax.plot(nu_cen, f3, '-', color='k', )

    def alpha_calc(self, nu, flux, errors):
        X = np.log10(np.array(nu))
        Y = np.log10(np.array(flux))
        # alpha = X[0]-
        sol = np.polyfit(X, Y, 1)
        alpha = sol[0]
        # B = sol[1]
        sigma_alpha = np.sqrt((errors[0] ** 2) * (
                1.0 / (np.log10(nu[0] / nu[1])) * 1.0 / (
                np.log(10) * flux[0])) ** 2 + (errors[1] ** 2) * (
                                      1.0 / (
                                  np.log10(nu[0] / nu[1])) * -1.0 / (
                                              np.log(10) * flux[
                                          1])) ** 2)
        return [alpha, sigma_alpha]

    def par_string(self, fit_par, perr, par_names, funcname):
        if fit_par!=None and perr !=None:
            mystring = '\n' + funcname + '\n'
            for par, err, name in zip(fit_par, perr, par_names):
                pstring = name +r'= {:.2f} $\pm$ {:.2f}'.format(par, err) +'\n'
                mystring += pstring
        else:
            mystring = '\n' + funcname + '\n None'
        return mystring

    def tick_offsets(self, f, ax_cord):
        xp = f._ax1.xaxis.get_majorticklocs()
        yp = f._ax1.yaxis.get_majorticklocs()
        tick_pos_pix = list(zip(xp, yp))
        tick_pos_world = []
        dra = []
        ddec = []
        for i in tick_pos_pix:
            dum = f.pixel2world(i[0], i[1])
            sky = SkyCoord(ra=dum[0], dec=dum[1], unit=(u.deg, u.deg),
                           frame='fk5')
            tick_pos_world.append(sky)
            dum_dra, dum_ddec = sky.spherical_offsets_to(ax_cord)
            dra.append(np.round(dum_dra.to(u.arcsec).value, 2))
            ddec.append(np.round(dum_ddec.to(u.arcsec).value, 2))
        return dra, ddec

    ##########################################################################

    def func_iau_names(self, spos):
        niau = 'J'
        sign = '+'
        if spos.find(sign) == -1:
            sign = '-'
        pos = spos.split(sign)
        niau += pos[0][:4]
        niau += sign
        niau += pos[1][:2]
        return niau

    def chisq(self, func, xdata, ydata, popt, sigma):
        res = ydata - func(xdata, *popt)
        return np.sum((res / sigma) ** 2)


    def prep_fit_arr(self, nu_arr, flux_arr, eflux_arr, alpha):
        if len(alpha) > 0 and alpha[1] > -9999:
            nu_arr = np.append(nu_arr, alpha[0])
            flux_arr = np.append(flux_arr, alpha[1])
            eflux_arr = np.append(eflux_arr, alpha[2])
        return nu_arr, flux_arr, eflux_arr

    def generate_model(self, nu_arr, *pars):
        mypars = pars
        nu_arr = np.array(nu_arr)
        model_lists = {'PL': self.power_law, 'CPL': self.CPL,
                       'EFFA': self.EFFA_func, 'IFFA': self.IFFA_func,
                       'SSA': self.SSA_func, 'GCV': self.gen_curve_func}
        if self.model_name not in model_lists.keys():
            raise Exception('Invalid Model: Please check the model code')

        model_func = model_lists[self.model_name]
        flux = np.copy(nu_arr)
        pnu_msk = tuple([nu_arr > 0])
        flux[pnu_msk] = model_func(nu_arr[pnu_msk], *mypars)
        nv = -1.0 * nu_arr[~pnu_msk[0]]
        n_low = nv * 0.99
        n_up = nv * 1.01
        flux[~pnu_msk[0]] = np.log10(
            model_func(n_up, *mypars) / model_func(n_low,
                                                   *mypars)) / np.log10(
            n_up / n_low)
        return flux

    def chi_sq_fit(self, nu_arr, flux_arr, eflux_arr, guess_pars, model_name):
        self.model_name = model_name
        model_lists = {'PL': self.power_law, 'CPL': self.CPL,
                       'EFFA': self.EFFA_func, 'IFFA': self.IFFA_func,
                       'SSA': self.SSA_func, 'GCV': self.gen_curve_func}

        if guess_pars == None or len(guess_pars) == 0:
            s0 = flux_arr[np.argwhere(nu_arr == 1.4)]  # getting NVSS
            alpha = -1.0
            nu_t = np.mean(nu_arr[nu_arr > 0])
            if self.model_name == 'PL':
                guess_pars = [s0, alpha]
            else:
                guess_pars = [s0, alpha, nu_t]
        try:
            self.fit_res, cov1 = curve_fit(self.generate_model,
                                                          nu_arr, flux_arr,
                                                           guess_pars,eflux_arr)
            self.perr_res = np.sqrt(np.diag(cov1))
            mymsk = tuple([nu_arr>0])
            chi_res = self.redchisq_calc(model_lists[self.model_name],
                                     nu_arr[mymsk], flux_arr[mymsk], 
                                     eflux_arr[mymsk], self.fit_res)           
        except:
            self.fit_res = [-999]*len(guess_pars)
            chi_res =-9999
            self.perr_res = [-999]*len(guess_pars)

        return self.fit_res, self.perr_res, chi_res

    def redchisq_calc(self, model, nu, flux, err, pars):
        guess_pars = pars
        dof = max((len(nu) - len(guess_pars) - 1), 1)
        red_chi = (1.0 / dof)* np.sum(np.square((model(nu, *guess_pars) - flux)/err))
        return red_chi

    def estimate_guess_pars(self, freq_arr, flux_arr, model):
        #calculate alpha_high and alpha_low:
        nuh = np.max(freq_arr)
        nul = np.min(freq_arr[freq_arr>0])
        
        sh = flux_arr[np.argmax(freq_arr)]
        sl = flux_arr[np.argmin(freq_arr[freq_arr>0])]
        #calculate median index
        ind = np.argwhere(freq_arr==1.4)[0][0]
        if ind >-1 and np.min(freq_arr)<1.3:
            sm  = flux_arr[ind]
            num = freq_arr[ind]
        else:
            ind = np.argwhere(freq_arr==np.median(freq_arr))[0][0]
            sm  = flux_arr[ind]
            num = freq_arr[ind]
        
        if model == 'PL':
           alpha, s0 = np.polyfit(freq_arr, flux_arr, 1)
           guess_pars = [s0, alpha]
        else:
            al_h = np.log10(sh/sm)/np.log10(nuh/num)
            al_l = np.log10(sl/sm)/np.log10(nul/num)
            s0  = sm 
            alpha = al_h
            nu_p = num
            q =  -alpha/(2*np.log(nu_p))
            if al_h>0 and al_l>0:
                s0 = sh
                nu_p = nuh
            if al_h>0 and al_l<0:
                alpha =al_l
                q =  -alpha/(2*np.log(nu_p))
            if al_h<0 and al_l<0:
                alpha = np.mean([al_h, al_l])
                q = 0.0
            guess_pars = [s0, alpha, nu_p]
            if model == 'CPL':
                guess_pars = [s0, alpha, q]
            if model == 'GCV':
                guess_pars = [s0, al_l, al_h, nu_p]
        return guess_pars
            