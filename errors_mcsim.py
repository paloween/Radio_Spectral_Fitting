#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 10:36:00 2018

@author: pallavipatil
"""

'''
This will be a place to write Monte Carlo simulations to estimate errors on 
parameters of a given function. All the functions will be defined in a separate 
module and this will be imported to various python functions and scripts. I will 
define a toy model function for testing this code. 

'''
import sys
sys.path.append('../Radio-SED-Fitting')
import numpy as np
from scipy.stats import skew, kurtosis
import Radio_Models_func 

rmfit = Radio_Models_func.RadioModelFit()

def test_func_sindex(flux):
    n = len(flux)
    if n%2 != 0 :
        alpha =  -9999
    else:
        X = np.array(np.log10(flux[int(n/2):]))
        Y = np.array(np.log10(flux[:int(n/2)]))
        sol = np.polyfit(X,Y,1) # Fit a spectral index
        alpha = sol[0]
    return alpha


def normal_Gauss(x, mu, sigma):
    xn = (x-mu)/sigma
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-xn*xn/2)

def erf_gauss(x, mu, sigma):
    erf = np.zeros(len(x))
    nm_Gauss = normal_Gauss(x, mu, sigma)
    for ii in range(1,len(x)):
         erf[ii] = nm_Gauss[ii]+ erf[ii-1]
    return erf
        
    
    
def errors_mcsim(func, nvars1, nvars2, var, verr, ntrial, summ_F=False, plotF=False):
    '''
    vars1: total number of variables
    vars2: total number of variables to vary
    var: variables with the ones to vary at the beginning
    '''
    nerrf = 1000  # Size of the error function    
    #Construct an error function for gauss with mean 500 and sigma 100
    xmean = 500.0
    sigma = 100.0
    gvars = np.zeros(len(var))
    carlo = np.zeros(ntrial)
    # Create a normalized Gaussian error function
    erfx = np.linspace(0, nerrf, num=nerrf)
    erf = erf_gauss(erfx, xmean, sigma )
    nm_erfx = (erfx-xmean)/sigma
    
    statsum = []
    ssumm_names = ['Central Value', 'Mean', 'Avg. Deviation', 'Std dev', 'Variance',
                   'Skew', 'Kurtosis', 'Median', 'e1Low', 'e1High', 'e2Low', 'e2High']
    # Optimum value from central values/measured/estimated values of parameters
    
    #gvars array will hold randomly selected variable values based on gaussian 
    #distribution around the central value. If the variable fixed then the fixed 
    # value. 
    opv = func(*var)
    statsum.append(opv)
    if nvars1 != nvars2:
        for jj in range(nvars2, nvars1):
            gvars[jj] = var[jj]
    
    for nt in range(ntrial):
        # means : var
        # sigs : verr
        randg = np.random.randn(nvars2)
        gvars[:nvars2] = var[:nvars2] + randg*verr[:nvars2]
        carlo[nt] = func(*gvars)
    
    mean = np.nanmean(carlo)
    adev = np.nansum(carlo-mean)/len(carlo)
    rms = np.nanstd(carlo)
    std2 = np.nanvar(carlo)
    skw = skew(carlo, nan_policy='omit')
    curt = kurtosis(carlo, nan_policy='omit')
    statsum.extend([mean, adev, rms, std2, skw, curt])
    quants = np.nanpercentile(carlo,[50, 15.87, 84.13, 2.27, 97.72] )
    statsum.extend(quants)
    if summ_F: return dict(zip(ssumm_names, statsum))
    if plotF: return carlo
    else:
        return statsum
    
    
   
def bootstrap_mcsim(xarr, yarr, yerr,  guess_pars, model,ntrial):             
    nerrf = 1000  # Size of the error function    
    #Construct an error function for gauss with mean 500 and sigma 100
    xmean = 500.0
    sigma = 100.0
    carlo = np.zeros((ntrial, len(guess_pars)))
    calculate_par = False
    # Create a normalized Gaussian error function
    '''
    erfx = np.linspace(0, nerrf, num=nerrf)
    erf = erf_gauss(erfx, xmean, sigma )
    nm_erfx = (erfx-xmean)/sigma
    '''
    par_bestfits = []
    ssumm_names = [ 'Mean', 'Avg. Deviation', 'Std dev', 'Variance',
                   'Skew', 'Kurtosis', 'Median', 'e1Low', 'e1High', 'e2Low', 'e2High']
    # Optimum value from central values/measured/estimated values of parameters
    
    #gvars array will hold randomly selected variable values based on gaussian 
    #distribution around the central value. If the variable fixed then the fixed 
    # value.
    '''
    opv = func(*var)
    statsum.append(opv)
    if nvars1 != nvars2:
        for jj in range(nvars2, nvars1):
            gvars[jj] = var[jj]'''
    nvars = len(yarr)
    var = yarr
    verr = yerr
    if model == 'CPL':
        calculate_par = True
        func = rmfit.CPL
    sp_carlo = []
    nup_carlo = []
    for nt in range(ntrial):
        # means : var
        # sigs : verr
        randg = np.random.randn(nvars)
        gvars = var+ randg*verr
        fit_res, perr_res, chi_res = rmfit.chi_sq_fit(xarr, gvars,
                         np.abs(gvars)*0.005, guess_pars, model)
        for ii in range(len(fit_res)):
            carlo[nt][ii]= fit_res[ii]
        if calculate_par:
            nupeak = np.exp(-fit_res[1]/(2*fit_res[2]))
            speak = func(nupeak, *fit_res)
            sp_carlo.append(speak)
            nup_carlo.append(nupeak)
        
    
    tcarlo = np.transpose(carlo)
    for kk in range(len(fit_res)):
        statsum = []
        pararr = tcarlo[kk][:]
        mean = np.mean(pararr)
        adev = np.sum(pararr-mean)/len(pararr)
        rms = np.std(pararr)
        std2 = np.var(pararr)
        skw = skew(pararr)
        curt = kurtosis(pararr)
        statsum.extend([mean, adev, rms, std2, skw, curt])
        quants = np.percentile(pararr,[50, 15.87, 84.13, 2.27, 97.72] )
        statsum.extend(quants)
        par_bestfits.append(statsum)
    if calculate_par:
        sp_carlo= np.array(sp_carlo)
        nup_carlo  = np.array(nup_carlo)
        for pararr in [sp_carlo, nup_carlo]:
            statsum = []
            mean = np.mean(pararr)
            adev = np.sum(pararr-mean)/len(pararr)
            rms = np.std(pararr)
            std2 = np.var(pararr)
            skw = skew(pararr)
            curt = kurtosis(pararr)
            statsum.extend([mean, adev, rms, std2, skw, curt])
            quants = np.percentile(pararr,[50, 15.87, 84.13, 2.27, 97.72] )
            statsum.extend(quants)
            par_bestfits.append(statsum)
                
    return ssumm_names, par_bestfits, tcarlo




def run_test():
    
    #errors_mcsim(func, nvars1, nvars2, var, verr, ntrial, summ_F=False):
    #test_func = test_func_sindex(flux)
    
    #flux and freq
    
    flux = [7.7678, 6.5996, 8.599000, 11.399]
    nvars1 = 4
    nvars2 = 2
    ntrial = 10000
    verr = [0.0471, 0.043]
    summF = True
    
    test_sum = errors_mcsim(test_func_sindex, nvars1, nvars2, flux, verr, ntrial, 
                            summF, False)
    
    for key, value in zip(test_sum.keys(), test_sum.values()):
        print ('{} : {:.5f} '.format(key, value ))
        
        
            
    
    
    
    
    





    