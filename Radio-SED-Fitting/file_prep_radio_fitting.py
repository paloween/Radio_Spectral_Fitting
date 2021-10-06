# -*- coding: utf-8 -*-
"""

@author: ppatil
"""
import numpy as np
# import os
from astropy.table import Table, join
import pandas as pd
import sys
sys.path.insert(1,'../Radio-SED-Fitting')
import config


class LoadTables:
    def __init__(self):
        self.atscat = Table.read(config.WISE_XMATCH_FILE)
        self.atgcat = Table.read(config.GLEAM_FILE)
        self.vla_ax = Table.read(config.JMFIT_AX_FILE, format='csv',
                            delimiter=',')
        self.vla_bx = Table.read(config.JMFIT_BX_FILE, format='csv',
                            delimiter=',')

        #self.vla_ax_grp = vla_ax.group_by(keys='WISEname')
        #self.vla_bx_grp = vla_bx.group_by(keys='WISEname')
        #self.vla_ax_grp.add_index('WISEname')
        #self.vla_bx_grp.add_index('WISEname')
        self.atscat.add_index('WISEname')
        self.atgcat.add_index('WISEname')

        sp_df = pd.read_excel(config.SP_CLASS_FILE)
        self.sp_class = Table.from_pandas(sp_df)

        snr_ax = Table.read(config.SNR_AX_FILE, format='csv')
        snr_bx = Table.read(config.SNR_BX_FILE, format='csv')

        self.sp_class = join(self.sp_class, snr_ax, join_type='left',
                             keys='Source_name')
        self.sp_class = join(self.sp_class, snr_bx, join_type='left',
                             keys='Source_name')
        
        #Change on Feb 10
        # Adding VCSS catalog
        self.vcss_cat = Table.read(config.VLITE_FILE)
        self.vcss_cat.add_index('Source_name')
    def get_tabs(self):
        return self.atscat, self.atgcat, self.vla_ax, self.vla_bx

    def get_sp_info(self):
        return self.sp_class


class FilePrepUtils:

    def __init__(self):

        self.arx_freqs = [1.4, 1.4, 0.843,  0.150, 0.074, 0.325, 0.352, 4.85,
                          0.338, 3.0]
        self.arx_fcols = ['FNVSS', 'FPEAK', 'Fsumss', 'Ftgss', 'Fvlssr',
                          'Speak_WENSS', 'Fwish', 'Fgb6', 'FVLITE_p', 'Total_flux']
        self.arx_ecols = ['ENVSS', 'RMS','Esumss', 'Etgss', 'Evlssr',
                          'Nse_WENSS', 'Ewish', 'Egb6','EVLITE_p', 'E_Total_flux']
        gleam_freq = [76, 84, 92, 99, 107, 115, 122, 130, 143, 151, 158, 166,
                      174, 181, 189, 197, 204, 212, 220, 227]  # MHz
        self.gleam_f_cols = []
        self.gleam_ferr_cols = []
        self.cat_names = ['NVSS', 'FIRST', 'SUMSS', 'TGSS', 'VLSSr', 'WENSS',
                          'WISH', 'GB6', 'VLITE', 'VLASS']
        self.oir_fcols = ['FWISE1', 'FWISE2', 'FWISE3', 'FWISE4', 'F70', 'F170',
                          'F250', 'F350', 'F500']
        self.oir_ecols = ['EWISE1', 'EWISE2', 'EWISE3', 'EWISE4', 'E70', 'E170',
                          'E250', 'E350', 'E500']
        self.oir_freq = [88174, 66620, 24982, 12491, 4282, 1763, 1199, 856, 599]

        for frq in gleam_freq:
            if frq < 100:
                self.gleam_f_cols.append('peak_flux_0' + str(frq))
                self.gleam_ferr_cols.append('err_peak_flux_0' + str(frq))
            else:
                self.gleam_f_cols.append('peak_flux_' + str(frq))
                self.gleam_ferr_cols.append('err_peak_flux_' + str(frq))

    '''
    Note: Added 3%flux error in this. The error on spectral index does include 
    the 3% systematic error. 
    If you change the errors in the FSB file. Remove this part. 
    Feb 8 2019: removing the 3% error
    '''

    def check_sourcestruct(self, row):
        if not np.ma.is_masked(row['P_flux']):
            if row['result'] == 'UR':
                flux = row['P_flux']
                flux_err = row['pflux_err']
                ferr = flux_err
            else:
                flux = row['I_flux']
                flux_err = row['iflux_err']
                ferr = flux_err
        else:
            flux = -999.0
            ferr = -999.0
        return flux, ferr

    def check_Spindex(self, row):
        if not np.ma.is_masked(row['SpIdx']):
            spid = row['SpIdx'][0]
            sperr = row['SpIdx_err'][0]
        else:
            spid = -9999.0
            sperr = -9999.0
        return spid, sperr

        # modifying Data prep to add VLITE data point

    def data_prep(self, wisen, glmcat, scat, jvla_AX, jvla_BX):
        path_tab = '/Users/ppatil/Desktop/VLA/Radio_Fits_v2/New_Analysis/'
        final_comptab = Table.read(path_tab+'final_tab_new_jvla.csv')
        final_comptab.add_index('Source_name')
        
        vcss_cat = Table.read(config.VLITE_FILE)
        vcss_cat.add_index('WISEname')
        racs_cat = Table.read(config.RACS_FILE)
        racs_cat.add_index('WISEname')
        # These array are for the radio fitting routine:
        flux_arr = []
        freq_arr = []
        eflux_arr = []

        freqs = []
        fluxes = []
        errors = []
        labels = []

        alpha_AX = []
        alpha_BX = []
        alpha_GL = []

        freq_OIR = []
        flux_OIR = []
        eflux_OIR = []

        FAX_10GHZ = 0
        EAX_10GHZ = 0
        FBX_10GHZ = 0
        EBX_10GHZ = 0

        rarx_flag = 0
        fgl = 0
        egl = 0

        arx_fluxes = []
        arx_errors = []

        '''
         = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisename]
         = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisename]   
        '''

        for reg in jvla_AX:
            [f, ferr] = self.check_sourcestruct(reg)
            if f > -999.0:
                FAX_10GHZ += f
                EAX_10GHZ += ferr

        for reg in jvla_BX:
            [f, ferr] = self.check_sourcestruct(reg)
            if f > -999.0:
                FBX_10GHZ += f
                EBX_10GHZ += ferr

        for fcol, ecol in zip(self.arx_fcols, self.arx_ecols):
            flux = scat[fcol]
            ferr = scat[ecol]
            if fcol == 'Ftexas':
                flux = flux * 1000
                ferr = ferr * 1000
            if fcol == 'Fvlssr':
                flux = flux * 1000
                ferr = ferr * 1000
            if fcol == 'Total_flux': # VLASS
                flux = flux*1000
                ferr = flux*0.2

            arx_fluxes.append(flux)
            arx_errors.append(ferr)

        for jj in range(len(arx_fluxes)):
            if arx_fluxes[jj] > 0:
                fluxes.append(arx_fluxes[jj])
                errors.append(arx_errors[jj])
                freqs.append(self.arx_freqs[jj])
                labels.append(self.cat_names[jj])

        ''' A modification : Date May 3rd 2018::
        REplacing all of the GLEAM observations with a single wide band flux 
        point and a spectral index'''

        # glmcat = atgcat.loc[wisename]
        if not np.isnan(glmcat['peak_flux_wide']):
            fgl = float(glmcat['peak_flux_wide'] * 1000)
            egl = float(glmcat['err_peak_flux_wide'] * 1000)
            fluxes.append(fgl)
            errors.append(egl)
            freqs.append(0.2005)
            labels.append('GLEAM')
        if not np.isnan(glmcat['alpha']):
            alpha_GL = [-0.2005, glmcat['alpha'], glmcat['err_alpha']]

        freqs_arr = np.array(freqs)
        fluxes_arr = np.array(fluxes)
        errors_arr = np.array(errors)

        # Other radio observations:
        rarx_flux = []
        rarx_ferr = []
        rarx_freq = []
        rad_frq = [32, 10.45, 8.4, 5, 4.9, 4.85, 4.77, 2.7, 1.4, 0.843, 0.408,
                   0.365, 0.325, 0.330, 0.325, 0.160, 0.151, 0.08, 0.074]
        # Look for index of the first freq

        mir_row = glmcat
        r_ind = mir_row.colnames.index('F32GHZ')
        re_ind = mir_row.colnames.index('E32')
        r_cols = mir_row.colnames[r_ind: r_ind + len(rad_frq) + 1]
        e_cols = mir_row.colnames[re_ind: re_ind + len(rad_frq) + 1]

        r_cols.remove('NVSS1.4')
        e_cols.remove('ENVSS1.4')

        FALMA = mir_row['F870']
        EALMA = mir_row['E870']

        for fcol, ecol, nu in zip(self.oir_fcols, self.oir_ecols,
                                  self.oir_freq):
            if mir_row[fcol] > 0:
                freq_OIR.append(nu)
                flux_OIR.append(mir_row[fcol])
                eflux_OIR.append(mir_row[ecol])
        OIR = [freq_OIR, flux_OIR, eflux_OIR]

        for fcol, ecol, nu in zip(r_cols, e_cols, rad_frq):
            if mir_row[fcol] > 0 and isinstance(mir_row[ecol], float):
                rarx_flux.append(mir_row[fcol])
                rarx_ferr.append(mir_row[ecol])
                rarx_freq.append(nu)

        if len(rarx_flux) > 0:
            rarx_flag = 1
        flux_arr.extend(fluxes_arr)
        freq_arr.extend(freqs_arr)
        eflux_arr.extend(errors_arr)
        # Other radio archival data
        if rarx_flag == 1:
            for nu, fnu, efnu in zip(rarx_freq, rarx_flux, rarx_ferr):
                if ~np.isin(nu, freq_arr):
                    freq_arr.append(nu)
                    flux_arr.append(fnu)
                    eflux_arr.append(efnu)
                    labels.append('RARX')

        # JVLA Observations:

        # Working on A observations:
        if FAX_10GHZ > 0:
            flux_arr.append(FAX_10GHZ)
            freq_arr.append(10)
            eflux_arr.append(EAX_10GHZ)
            labels.append('AX')
        if FAX_10GHZ > 0 and len(jvla_AX) == 1:
            [spidx, spidx_e] = self.check_Spindex(jvla_AX)
            nu_0 = 10
            if spidx>-999:
                source = jvla_AX['Source_name']
                finalt = final_comptab.loc[source]
                if len(finalt)<5:
                    finalt = finalt[finalt['regionAX']=='reg1']
                spidx = finalt['alIB_AX']
                spidx_e = finalt['alIB_err_AX']
            alpha_AX = [-nu_0, spidx, spidx_e]
        if FBX_10GHZ > 0:
            flux_arr.append(FBX_10GHZ)
            freq_arr.append(10)
            eflux_arr.append(EBX_10GHZ)
            labels.append('BX')
        if FBX_10GHZ > 0 and len(jvla_BX) == 1:
            [spidx, spidx_e] = self.check_Spindex(jvla_BX)
            nu_0 = 10
            if spidx>-999:
                source = jvla_BX['Source_name']
                finalt = final_comptab.loc[source]
                if len(finalt)<5:
                    finalt = finalt[finalt['regionBX']=='reg1']
                spidx = finalt['alIB_BX']
                spidx_e = finalt['alIB_err_BX']
            alpha_BX = [-nu_0, spidx, spidx_e]

        flux_arr = np.array(flux_arr)
        freq_arr = np.array(freq_arr).astype('float')
        eflux_arr = np.array(eflux_arr).astype('float')
        int_arr = []
        for ii, nu in enumerate(freq_arr):
            fl = flux_arr[ii]
            efl = eflux_arr[ii]
            if fl <= 0 or efl <= 0 or np.isnan(efl):
                int_arr.append(ii)

        freq_arr = np.delete(freq_arr, int_arr)
        flux_arr = np.delete(flux_arr, int_arr)
        eflux_arr = np.delete(eflux_arr, int_arr)
        sp_flux = [FAX_10GHZ, EAX_10GHZ, FBX_10GHZ, EBX_10GHZ, fgl, egl]
    
        '''
        Change made on Feb 10 2021: Adding VCSS Catalog)
        '''
        #Remove the Older VLITE data and add new one from VCSS
        if 'VLITE' in labels:
            myind = labels.index('VLITE')
            freq_arr = np.delete(freq_arr, myind)
            flux_arr = np.delete(flux_arr, myind)
            eflux_arr = np.delete(eflux_arr, myind)
            labels = np.delete(labels, myind)
        
        
        try:
            vcss_row = vcss_cat.loc[wisen]
            #cmpt = vcss_row['compactness']
            #if cmpt<0.8:
            #    col = 'total_flux_1'
            #else: col = 'peak_flux_1'
            col = 'total_flux'
            vcss_flux = vcss_row[col]
            vcss_eflux = vcss_row['e_'+col]
            flux_arr = np.append(flux_arr, vcss_flux)
            eflux_arr = np.append(eflux_arr, vcss_eflux)
            freq_arr = np.append(freq_arr, 0.338)
            labels = np.append(labels, 'VLITE')
        except:
            pass
        
        #Read RACS info
        try:
            racs_row = racs_cat.loc[wisen]
            #cmpt = vcss_row['compactness']
            #if cmpt<0.8:
            #    col = 'total_flux_1'
            #else: col = 'peak_flux_1'
            col = 'total_flux_source'
            racs_flux = racs_row[col]
            racs_eflux = racs_row['e_'+col]
            flux_arr = np.append(flux_arr, racs_flux)
            eflux_arr = np.append(eflux_arr, racs_eflux)
            freq_arr = np.append(freq_arr, 0.887)
            labels = np.append(labels, 'RACS')
        except:
            pass
            
        
        return freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, [
            FALMA, EALMA], OIR, sp_flux, list(labels)
