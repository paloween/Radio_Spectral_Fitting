#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 20:25:04 2019

@author: pallavipatil
"""



import numpy as np
from astropy.cosmology import Planck15 as cm 
from scipy import polyfit
# calculate distances::
nuLow = 0.01
nuUpp = 100


def ang_diameter_dist(z):
    return cm.comoving_distance(z).value/(1+z)

def luminosity_dist(z):
    return cm.comoving_distance(z).value*(1+z)



def X_q_alpha (q, alpha):
    nu1 = nuLow
    nu2 = nuUpp
    return (nu2**(q-alpha) - nu1**(q-alpha))/(q-alpha)


def total_radio_luminosity(z,S_10GHz,  alpha_thin, nu):
    num1 = 1.2 * 10**(36) * (1+z)**(1+alpha_thin)*S_10GHz
    num2 = cm.comoving_distance(z).value**2 *(nu)**alpha_thin
    num3 = X_q_alpha(1, alpha_thin)
    L_rad = num1 * num2 * num3
    return L_rad # in erg/s

def spectral_radio_luminosity_old(z,nu,  S_nu, alpha, nu2):
    d_L = luminosity_dist(z)
    #nu2o = nu2/(1+z)
    Snu2 = S_nu*(nu2/nu)**(alpha)
    L_nu_r = 4*np.pi*d_L*d_L*Snu2*(1+z)**(-1-alpha)
    return L_nu_r


def spectral_radio_luminosity(z, nu, S_nu, alpha):
    d_L = luminosity_dist(z)
    L_nu = 4*np.pi*d_L*d_L*S_nu/((1+z)**(1+alpha))
    return L_nu


def min_mag_field(S_10GHz, z, alpha_thin, nu, ff_rl, theta_maj, theta_min,
                  a):
    num1 = S_10GHz/(theta_maj*theta_min)
    num2 = a * (1+z)**(4+alpha_thin) * nu**(alpha_thin)
    num3 = X_q_alpha(0.5, alpha_thin)
    num4 = ff_rl*theta_min*cm.comoving_distance(z).value
    Bmin = 2.93 *10**(-4)*(num1*num2*num3/num4)**(2/7)
    return Bmin #in gauss
    
def relativ_pressure(Bmin):
    return 0.031*Bmin*Bmin



def ssa_mag_field(nu_p, S_p, z, theta_maj, theta_min):
    # unit conversions:
    sp = S_p/1000.0 # In Jy
    thmaj = theta_maj*1000 #In mas
    thmin = theta_min*1000
    num = nu_p**5*(thmaj*thmin)**2
    dem = 8**5*sp**2*(1+z)
    H = num/dem
    return H


def total_lobe_energy(z,ff_rl, theta_maj, theta_min,Bmin ):

    num1 = theta_maj*theta_min*theta_min*ff_rl**(3/7)
    num2 = cm.comoving_distance(z)**3 * (1+z)**(-3)
    E_rel = 1.6*10**56*num1*num2*Bmin*Bmin  #Bmin is calculated at ff_rl = 1
    return E_rel

  
def t_half(Bmin, nu):
    return 0.68*Bmin**(-3/2)*nu**(-0.5) # in years

def electron_density(Bmin, a, alpha_thin):
    num1 = 4.20 * 10**3/a
    num2 = X_q_alpha(0,alpha_thin)/X_q_alpha(0.5, alpha_thin)
    ne = num1*num2*Bmin**(5/2)
    return ne

def ratio_rme_par_en(Bmin, a, alpha_thin):
    return 68/a * X_q_alpha(0,alpha_thin)/X_q_alpha(0.5, alpha_thin)* Bmin**(0.5)

def power_law_est(nu1, S1, nu2, alpha):
    return S1*(nu2/nu1)**(-alpha)

def two_pt_alpha(freq, flux, fluxerr, s0_flag=False):
    x = freq
    y = flux
    errors = fluxerr
    X = np.log10(np.array(x))
    Y = np.log10(np.array(y))
    sol = polyfit(X,Y,1)
    alpha = sol[0]
    sigma_alpha = np.sqrt((errors[0]**2)*(1.0/(np.log10(x[0]/x[1]))*1.0/(np.log(10)*y[0]))**2 + (errors[1]**2)*(1.0/(np.log10(x[0]/x[1]))*-1.0/(np.log(10)*y[1]))**2)
    if s0_flag:
        return alpha, sigma_alpha, 10**sol[1]
    else:
        return alpha, sigma_alpha


def two_pt_alpha_mc(*flux):
    n = len(flux)
    if n%2 != 0 :
        alpha =  -9999
    else:
        X = np.array(np.log10(flux[int(n/2):]))
        Y = np.array(np.log10(flux[:int(n/2)]))
        sol = polyfit(X,Y,1) # Fit a spectral index
        alpha = sol[0]
    return alpha

