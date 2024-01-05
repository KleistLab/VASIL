#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 16:28:56 2023

@author: raharinirina
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:53:58 2023

@author: raharinirina
"""
import numpy as np
import pandas as pd
import numpy.ma as ma
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import pdb
import seaborn as sns
import pickle
import mpl_axes_aligner
import re

"""
Requires system arguments

sys.argv[1]: "path/to/Country_1/new/path/to/Country_2/new/path/to/Country_3/..." ### as many countries as desired but always separated by /new/ and folder must include all country results
sys.argv[2]: "Country_1/Country_2/Country_3/..." ### country names separated by / 
sys.argv[2]: "Country_1/Country_2/Country_3/..." ### country labels separated by /
sys.argv[4]: percentage to filter spikegroups #e.g. 0.01 
sys.argv[5]: dates to start the comparison for each variants, separated by /, if only one then they will be all the same ### e.g. "2022-04-09/2022-04-26" for 2 variants
sys.argv[6]: dates to end the comparison for each variants, separated by /, if only one then they will be all the same ### e.g. "2022-08-21/2023-04-15" for 2 variants
sys.argv[7]: integer number of variants to check
sys.argv[8]: folder path to store the results 
sys.argv[9:9+sys.argv[7]]: list all lineages to simulate ### e.g "BE.10.ALL", "BE.10.ALL" 
sys.argv[9+sys.argv[7]:len(sys.argv)]: list 1 color for each countries by their order in sys.argv[2] ### e.g. "red" "orange" for two countries
"""

def moving_average(X, window = 7):
    u = np.zeros(len(X))
    u[:window] = X[:window]
    for i in range(window, len(X)):
        u[i] = np.mean(X[i-window:i+1])
    
    return u


import matplotlib
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)

   
def plot_fit(countries_dir_list, countries_list, countries_labels, lineage_list, color_list, w_save = len(sys.argv)-1):
    status_list = []
    status_list_pseudo =[]
    lineage_list_countries = []    
    for k in range(len(lineage_list)):
        t_dates_list = []
        S_all_mean_list = []
        all_dates = []
        all_prop_dates = []
        lineage_freqs_list = []
        prop_mask_list = []
        day_prop_list = []
        for c_ind in range(len(countries_list)):
            S_mean_file = countries_dir_list[c_ind]+"/results/Susceptible_weighted_mean_over_spikegroups_all_PK.csv"
            # needs to be updated to allow individual weighting 
            S_mean_df = pd.read_csv(S_mean_file)
            S_all_mean_list.append(S_mean_df.to_numpy()[:, S_mean_df.columns != "Days"].astype(float))
            t_dates_list.append(S_mean_df["Days"].tolist())
            all_dates += list(S_mean_df["Days"])
            
            # processing of lineage frequency data            
            file = open(countries_dir_list[c_ind]+"/Spikegroups_membership.pck", "rb") #hard-coded for vaccination simulations
            Pseudogroup_dic = pickle.load(file)
            file.close() 
            
            lineage_freq = Pseudogroup_dic["Frequencies"].copy()
            threshold = float(sys.argv[4])
            try:
                lineage_freq.drop(columns = "Unnamed: 0", inplace = True)
            except:
                pass
        
            ### Drop proportions column that does not start with correspond to E[Susceptible] files
            lineage_freq = lineage_freq.drop(index = lineage_freq.index[:list(lineage_freq['date']).index(S_mean_df["Days"][0])])
            freqs = lineage_freq.loc[:, lineage_freq.columns != 'date']
            # imputing frequencies below threshold and normalization
            prop_mask = np.all(lineage_freq.loc[:, lineage_freq.columns != 'date'] == 0.0, axis = 1)
            freqs = freqs.mask(freqs < threshold)
            freqs = freqs.fillna(0)
            col_sums = freqs.sum(axis = 1).values
            freqs = freqs.divide(col_sums, axis="rows")
            freqs = freqs.fillna(0)
            lineage_freq.loc[:, lineage_freq.columns != 'date'] = freqs
            day_prop = lineage_freq["date"].tolist()
            
            all_prop_dates += day_prop
            day_prop_list.append(day_prop)
            lineage_freqs_list.append(lineage_freq)
            prop_mask_list.append(prop_mask)
            
        
        # plotting
        PreFig(xsize = 20, ysize = 20)
        fig = plt.figure(figsize = (15, 7))
        ax = fig.add_subplot(1, 1, 1)
        ### end of observation line
        all_dates = list(np.unique(all_dates))
        all_dates.sort(key = lambda date: datetime.strptime(date, "%Y-%m-%d")) 
        
        all_prop_dates = list(np.unique(all_prop_dates))
        all_prop_dates.sort(key = lambda date: datetime.strptime(date, "%Y-%m-%d")) 
        
        if all_prop_dates[-1]>all_dates[-1]:
            indx = all_prop_dates.index(all_prop_dates[-1])
            all_dates = all_dates + [all_prop_dates[indx+i] for i in range(len(all_prop_dates[indx:]))]
            
        t_prop_all = np.arange(len(all_prop_dates)).astype(int)
        inds_dates_all = np.arange(len(all_dates)).astype(int)
        
        #ax.axvline(x = len(all_dates) - 1, ymin = -1, ymax = 1, ls = "--", linewidth = 2, color = "grey")
        # different axis for proportions
        ax_twin = ax.twinx()
        
        ### Separate figure for relative fitness vs change in proportion
        PreFig(xsize = 20, ysize = 20)
        fig2 = plt.figure(figsize = (15, 7))
        ax2 = fig2.add_subplot(1, 1, 1)
        # different axis for proportions
        ax2_twin = ax2.twinx()
        
        ### Ticks are the same everywhere
        try:
            start = str(sys.argv[5]).split("/")
            if len(start) == 1:
                start = start[0]
            else:
                start = start[k]
            
            end = str(sys.argv[6]).split("/")
            if len(end) == 1:
                end = end[0]
            else:
                end = end[k]
            
            if start in list(all_dates):
                x_min = list(all_dates).index(start)
            else:
                wd = 0
                notstop = True
                while wd < len(all_dates) and notstop:
                    if all_dates[wd]<start:
                        notstop = True
                        wd +=1
                    else:
                        notstop = False    
                x_min = wd
                
            if end not in list(all_dates):
                x_max = (len(all_dates) - 1) 
            else:
                x_max = list(all_dates).index(end)
        except:
            x_min = None
    
        if (x_min is not None):                
            t_dates_show = np.array(all_dates)[x_min:x_max+1]
            check_last = x_max
            if len(t_dates_show)>200:
                pp = 7*4
            else:
                pp = min(len(t_dates_show), 14)
            perday = inds_dates_all[x_min:x_max+1][::pp]

        else:
            t_dates_show = all_dates
            check_last = len(all_dates) - 1
            if len(t_dates_show)>200:
                pp = 7*4
            else:
                pp = min(len(t_dates_show), 14)
                
            perday = inds_dates_all[::pp]
            
        date_ticks = np.array(all_dates)[perday].tolist()
        if all_dates[check_last] not in date_ticks:
            try:
                n=list(all_dates).index(date_ticks[-1])+pp
            except:
                n= perday[-1]+pp
            while n<len(all_dates)-1:
                date_ticks.append(all_dates[n])
                perday = np.append(perday, n)
                n += pp
                
            if check_last-perday[-1]<np.ceil(pp/5):
                perday = perday[:-1]
                date_ticks = date_ticks[:-1]
                
            date_ticks.append(all_dates[check_last])
            perday = np.append(perday, check_last)  
        
        
        
        save_plot = []
        other_trends = []
        other_labs = []
        other_cols = []
        other_props = []
        other_t_all = []
        other_prop_all = []
        other_t = []
        for c_ind in range(len(countries_list)):
            
            ES_df_dir = countries_dir_list[c_ind]+"/results/Immunological_Landscape_ALL"
            day_prop = day_prop_list[c_ind]
            t_prop = t_prop_all[list(all_prop_dates).index(day_prop[0]):list(all_prop_dates).index(day_prop[-1])+1]
            prop_mask = prop_mask_list[c_ind]
            t_dates = t_dates_list[c_ind]
            inds_dates = inds_dates_all[list(all_dates).index(t_dates[0]):list(all_dates).index(t_dates[-1])+1]
            S_all_mean = S_all_mean_list[c_ind]
            variant_freq = lineage_freqs_list[c_ind].copy()
            Pseudo_done_global = []
            
            file = open(countries_dir_list[c_ind]+"/Spikegroups_membership.pck", "rb") #hard-coded for vaccination simulations
            Pseudogroup_dic = pickle.load(file)
            file.close() 
            
            if c_ind == 0: 
                try:
                    Pseudogroup_dic.pop("Frequencies")
                except:
                    pass
                
                """### The first country is always the one we want to illustrate"""
                if lineage_list[k][-4:] == ".ALL" and not re.search("/", lineage_list[k]):
                    splited_var = []
                    for x in list(Pseudogroup_dic.keys()):
                        if (x[:len(lineage_list[k][:-4])] == lineage_list[k][:-4]):
                            if len(x) == len(lineage_list[k][:-4]):
                                splited_var.append(x)
                            elif len(x)>len(lineage_list[k][:-4]):
                                if x[len(lineage_list[k][:-4])] == ".":
                                    splited_var.append(x)
                            
                else:
                    splited_var = np.array(lineage_list[k].split("/"))
                    splited_var = splited_var[~(splited_var == "")]
                    splited_var = splited_var[~(splited_var == " ")]
             
                splited_var0 = np.array(splited_var).copy()
                splited_var = []
                for var in splited_var0:
                    if var[-4:] == ".ALL":
                        for x in list(Pseudogroup_dic.keys()):
                            if x[:len(var[:-4])] == var[:-4]:
                                if len(x) == len(var[:-4]):
                                    splited_var.append(x)
                                elif len(x)>len(var[:-4]):
                                    if x[len(var[:-4])] == ".":
                                        splited_var.append(x)
                    else:
                        splited_var.append(var)
                
                
                """### Get all the lineages belonging to the spikegroups relevant to lineage_list[k]"""
                relevant_groups = []
                for var in splited_var:
                    for x in list(Pseudogroup_dic.keys()):
                        if var in list(Pseudogroup_dic.keys()):
                            if Pseudogroup_dic[x] == Pseudogroup_dic[var]:
                                if x not in relevant_groups:
                                    relevant_groups.append(x)   
                splited_var = relevant_groups
            else:
                """### The group we compare is always the same as the list of lineages that are relevant to lineage_list[k] for the first country"""
                splited_var = []
                for var in relevant_groups:
                    if var in list(Pseudogroup_dic.keys()):
                        splited_var.append(var)
            
            num_avail = 0
            lab_done = {}
            lab_done[lineage_list[k]] = ""
            Pseudo_done = {}
            Pseudo_done[lineage_list[k]] = ""
            plot_prop = []
            Pseudo_Prop = np.zeros((len(t_prop)))
            start = 0
            ES_list = []
            prop_list = []
            for x in range(len(splited_var)):
                lineage = splited_var[x]
                try:
                    if lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                        ES_df = pd.read_csv(ES_df_dir+"/Susceptible_SpikeGroup_%s_all_PK.csv"%Pseudogroup_dic[lineage])
                        #print("-------------------------------------------------------------")
                        #print(countries_labels[c_ind], ES_df_dir+"/Susceptible_SpikeGroup_%s_all_PK.csv"%Pseudogroup_dic[lineage])
                        num_avail +=1 
                        
                        try:
                            ES_df.drop(columns = "Unnamed: 0", inplace = True)
                        except:
                            pass
                        
                        es_cols = ES_df.columns
                        ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                        ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))    
                        run = True
                except:
                    try:
                        if "Placeholder"+ lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                            ES_df = pd.read_csv(ES_df_dir+"/Susceptible_SpikeGroup_%s_all_PK.csv"%lineage)
                            num_avail +=1
                            try:
                                ES_df.drop(columns = "Unnamed: 0", inplace = True)
                            except:
                                pass
                            
                            es_cols = ES_df.columns
                            ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                            ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))
                            run = True
                    except:
                        try:
                            if lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                                ES_df = pd.read_csv("results/Immunological_Landscape_ALL/Susceptible_SpikeGroup_%s_all_PK.csv"%Pseudogroup_dic[lineage])
                                num_avail +=1
                                try:
                                    ES_df.drop(columns = "Unnamed: 0", inplace = True)
                                except:
                                    pass
                                
                                es_cols = ES_df.columns
                                ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                                ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))
                                run = True
                        except:
                            try:
                                if lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                                    ES_df = pd.read_csv("results/Immunological_Landscape/Susceptible_SpikeGroup_%s_all_PK.csv"%Pseudogroup_dic[lineage])
                                    num_avail +=1
                                
                                    try:
                                        ES_df.drop(columns = "Unnamed: 0", inplace = True)
                                    except:
                                        pass
                                    
                                    es_cols = ES_df.columns
                                    ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                                    ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))
                                    run = True
                            except:
                                 try:
                                     if "Placeholder"+ lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                                         ES_df = pd.read_csv("results/Immunological_Landscape/Susceptible_SpikeGroup_%s_all_PK.csv"%lineage)
                                         num_avail +=1
                                         try:
                                             ES_df.drop(columns = "Unnamed: 0", inplace = True)
                                         except:
                                             pass
                                            
                                         es_cols = ES_df.columns
                                         ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                                         ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))
                                        
                                         run = True
                                 except:
                                    print("Computation needed: Expected Susceptible file is not available for %s"%lineage, ES_df_dir)
                                    run = False
                                    
                # processing of Proportions data
                if run:   
                    if lineage in list(Pseudogroup_dic.keys()):
                        #lab_k = lineage + "*"+"/" 
                        lab_k = lineage + "/"
                        plot_prop.append(True)
                        if lineage in variant_freq.columns.astype(str):
                            #if Pseudogroup_dic[lineage] not in (Pseudo_done[lineage_list[k]].split("/")):
                            if lineage not in Pseudo_done[lineage_list[k]].split("/"):
                                sub_prop_lin = moving_average(variant_freq[lineage], window = 14)                                  
                                prop_list.append(sub_prop_lin)
                                Pseudo_Prop += sub_prop_lin
                                    
                                lab_done[lineage_list[k]] = lab_done[lineage_list[k]][:-1] + " + "+lab_k
                                if Pseudogroup_dic[lineage] not in Pseudo_done_global: ### Repeating Pseudogroups in variants combinations cannot be accounted twice in final proportion plot
                                    Pseudo_done_global.append(lineage)#(Pseudogroup_dic[lineage]) 
                            else:
                                lab_done[lineage_list[k]] += lab_k
                               
                        if Pseudogroup_dic[lineage] not in (Pseudo_done[lineage_list[k]].split("/")):
                            Pseudo_done[lineage_list[k]] += lineage + "/" #Pseudogroup_dic[lineage] + "/"
                        
                    else: 
                        lab_k = lineage+"/"
                        if lineage in variant_freq.columns.astype(str):
                            plot_prop.append(True)
                            if "Placeholder"+lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                                sub_prop_lin = moving_average(variant_freq[lineage], window = 14)
                                Pseudo_Prop += sub_prop_lin
                                prop_list.append(sub_prop_lin)

                                lab_done[lineage_list[k]] = lab_done[lineage_list[k]][:-1] + " + "+lab_k
                                if Pseudogroup_dic[lineage] not in Pseudo_done_global: ### Repeating Pseudogroups in variants combinations cannot be accounted twice in final proportion plot
                                    Pseudo_done_global.append(lineage)
                            else:
                                lab_done[lineage_list[k]] += lab_k
     
                        else:
                            plot_prop.append(False)
                            Pseudo_Prop += ma.masked_array(np.zeros(len(t_prop)), mask = np.ones(len(t_prop), dtype = bool))
                            if "Placeholder"+lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                                lab_done[lineage_list[k]] += " + "+lab_k
                                sub_prop_lin = ma.masked_array(np.zeros(len(t_prop)), mask = np.ones(len(t_prop), dtype = bool))
                                    
                                prop_list.append(sub_prop_lin)
                            else:
                                lab_done[lineage_list[k]] += " + "+lab_k
                        
                        if "Placeholder"+lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                            Pseudo_done[lineage_list[k]]+="Placeholder"+ lineage +"/"
                     
                    start +=1
                   
            lab_status = lab_done[lineage_list[k]][:-1]
            if not re.search(".ALL",lineage_list[k]):
                lab_k = lab_done[lineage_list[k]][:-1]
            else:
                lab_k = lineage_list[k]
            
            if lab_k[:3] == " + ":
                lab_k = lab_k[3:]
                
            if lab_status[:3] == " + ":
                lab_status = lab_status[3:]
            
            lab_status_regrouped = []
            lab_status_splited = lab_status.split(" + ")
            grouped = []
            lead = []
            for lb in range(len(lab_status_splited)):
                if lab_status_splited[lb] not in grouped:
                    for x in list(Pseudogroup_dic.keys()):
                        if lab_status_splited[lb] in list(Pseudogroup_dic.keys()):
                            if (x in lab_status_splited) and (Pseudogroup_dic[x] == Pseudogroup_dic[lab_status_splited[lb]]):
                                if Pseudogroup_dic[lab_status_splited[lb]] not in lead:
                                    if (x != lab_status_splited[lb]):
                                        lab_status_regrouped.append(lab_status_splited[lb]+"/"+x)
                                        grouped += [lab_status_splited[lb], x]
                                    else:
                                        lab_status_regrouped.append(x)
                                        grouped += [x]
                                        
                                    lead.append(Pseudogroup_dic[Pseudogroup_dic[lab_status_splited[lb]]])
                                else:
                                    indx = lead.index(Pseudogroup_dic[lab_status_splited[lb]])
                                    lab_status_regrouped[indx] += "/"+x
                                    grouped += [x]

            lab_status = "%s : %s"%(countries_labels[c_ind], lab_status)
            lab_status_pseudo = "%s: %s"%(countries_labels[c_ind], (" + ").join(lab_status_regrouped))
            if (lab_k != "" and num_avail != 0):
                if (num_avail == len(ES_list)):
                    #ES_ranges= np.mean(np.array(ES_list), axis = 0) # compute the mean
                    Props = np.array(prop_list)[:, :ES_list[0].shape[0]]
                    Props_normed = np.divide(Props, np.sum(Props, axis = 0)[np.newaxis, :], out = np.zeros(Props.shape), where = np.sum(Props, axis = 0)[np.newaxis, :]!=0)
                    ES_ranges = np.sum(np.array(ES_list)*Props_normed[:, :, np.newaxis], axis = 0) ## weighted mean
                    ES_ranges[np.sum(Props, axis = 0) == 0, :] = np.mean(np.array(ES_list), axis = 0)[np.sum(Props, axis = 0) == 0, :]
                else:
                    sys.exit("Loaded E[Susceptible] files were more than what is available for groups %s, recheck the loading process script/plotting/relative_advantage_groups.py Line 122-196"%lineage_list[k])
                
                # calculation of change in relative frequency from model
                Pseudo_Prop2 = Pseudo_Prop.copy()
                Pseudo_Prop2[Pseudo_Prop2 < 0.05] = 0
                Pseudo_Prop2 = list(Pseudo_Prop2)
                Pseudo_Prop_aligned = np.zeros(len(t_dates))

                day_prop_aligned = []
                t_prop_aligned = np.zeros(len(t_dates))
                prop_mask_aligned = np.zeros(len(t_dates)).astype(bool)
                gamma_prop = np.zeros(len(t_dates))
                SI_mask = np.zeros(len(t_dates)).astype(bool)
                for l in range(len(t_dates)):
                    #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
                    ## Dates are already set to always be to be successive
                    if t_dates[l] in day_prop:
                        w_l = list(day_prop).index(t_dates[l])
                        Pseudo_Prop_aligned[l] = Pseudo_Prop[w_l]
                        prop_mask_aligned[l] = prop_mask.tolist()[w_l]
                        day_prop_aligned.append(t_dates[l])
                        t_prop_aligned[l] = inds_dates[w_l]
                        try:
                            if Pseudo_Prop2[w_l] == 0 or Pseudo_Prop2[w_l+1] == 0:
                                gamma_prop[l] = float('nan')
                                SI_mask[l] == True
                            else:
                                gamma_prop[l] = (Pseudo_Prop2[w_l+1]/Pseudo_Prop2[w_l]) - 1
                        except:
                            gamma_prop[l] = float('nan')
                            SI_mask[l] = True
                    else:
                        gamma_prop[l] = float('nan')
                        SI_mask[l] = True
                
                ### Making sure to aligne proportions and incidence timelines
                if len(day_prop_aligned) < len(day_prop):
                    miss_days = [day for day in day_prop if day not in day_prop_aligned]
                    day_prop_aligned = list(day_prop_aligned) + [day_prop[list(day_prop).index(miss_days[i])] for i in range(len(miss_days))]
                    prop_mask_aligned = list(prop_mask_aligned) + [prop_mask.tolist()[list(day_prop).index(miss_days[i])] for i in range(len(miss_days))]
                    t_prop_aligned = list(t_prop_aligned) + [t_prop_all[list(all_prop_dates).index(miss_days[i])] for i in range(len(miss_days))]
                    Pseudo_Prop_aligned = list(Pseudo_Prop_aligned) + [Pseudo_Prop[list(day_prop).index(miss_days[i])] for i in range(len(miss_days))]
                    
                    day_prop_aligned = list(day_prop_aligned)
                    prop_mask_aligned = np.array(prop_mask_aligned)
                    t_prop_aligned = np.array(t_prop_aligned)
                    Pseudo_Prop_aligned = np.array(Pseudo_Prop_aligned)
                    
                # calculation of relative fitness
                gamma_SI = np.zeros((len(t_dates), ES_ranges.shape[1]))
                for i in range(ES_ranges.shape[1]):
                    S_x = ES_ranges[:, i]
                    S_mean = S_all_mean[:, i]
                    gamma_SI[:, i] = np.divide(S_x - S_mean, S_mean, out = S_x, where = (S_mean != 0)&(S_x!=0))
                
                # get min max gamma over PK at each timepoints
                SI_mask = np.array(SI_mask) + prop_mask_aligned[:len(inds_dates)] ### props are already aligned with indicence date
                gamma_SI_min, gamma_SI_max = np.min(gamma_SI, axis = 1), np.max(gamma_SI, axis = 1)
                gamma_SI_max = ma.masked_array(gamma_SI_max, mask = SI_mask)
                gamma_SI_min = ma.masked_array(gamma_SI_min, mask = SI_mask)
                Pseudo_Prop_masked = ma.masked_array(Pseudo_Prop_aligned, mask = prop_mask_aligned)
                
                if x_min is not None:
                    gamma_SI_min = gamma_SI_min[(inds_dates>=x_min)&(inds_dates<=x_max)]
                    gamma_SI_max = gamma_SI_max[(inds_dates>=x_min)&(inds_dates<=x_max)]
                    inds_dates_gamma = inds_dates[(inds_dates>=x_min)&(inds_dates<=x_max)]
                else:
                    inds_dates_gamma = inds_dates.copy()
                    
                if c_ind == 0:
                    ### only plot the fitness of the main country investigated
                    ax.set_title(lineage_list[k].replace("/", " + ") +" %s"%countries_labels[0][:3].upper(), fontsize = 32)
                    ax_twin.fill_between(inds_dates_gamma, gamma_SI_min, gamma_SI_max, color = color_list[c_ind], alpha = 0.3, label = "%s"%countries_labels[c_ind])
                    ax.plot(t_prop_aligned, 100*Pseudo_Prop_masked, linewidth = 6, color = color_list[c_ind] , label = "%s"%countries_labels[c_ind])
                    
                    inds_dates_gamma_0, gamma_SI_min_0, gamma_SI_max_0 = inds_dates_gamma.copy(), gamma_SI_min.copy(), gamma_SI_max.copy()
                    t_prop_aligned_0, Pseudo_Prop_masked_0 = t_prop_aligned.copy(), Pseudo_Prop_masked.copy()
                    
                    lim_0, lim_1 = ax_twin.get_ylim()
                    if lim_1 < 0:
                        lim_1 = 0.1*np.abs(lim_1)
                    if lim_0 > 0:
                        lim_0 = -0.1*np.abs(lim_0)
                    
                    if x_min is not None:
                        Prop_ref = Pseudo_Prop_aligned[:len(inds_dates)][(inds_dates>=x_min)&(inds_dates<=x_max)]
                    else:
                        Prop_ref = Pseudo_Prop_aligned.copy()[:len(inds_dates)]
                else:
                    other_trends.append((gamma_SI_min, gamma_SI_max))
                    other_labs.append(countries_labels[c_ind])
                    other_cols.append(color_list[c_ind])
                    if x_min is not None:
                        other_props.append(Pseudo_Prop_aligned[:len(inds_dates)][(inds_dates>=x_min)&(inds_dates<=x_max)])
                        other_t.append(inds_dates_gamma)
                    else:
                        other_props.append(Pseudo_Prop_aligned[:len(inds_dates)])
                        other_t.append(inds_dates_gamma)
                    
                    other_prop_all.append(Pseudo_Prop_aligned)
                    other_t_all.append(t_prop_aligned)
                        
                    #ax_twin.fill_between(inds_dates_gamma, gamma_SI_min, gamma_SI_max, color = color_list[c_ind], alpha = 0., label = "%s"%countries_labels[c_ind])
                    #ax.plot(t_prop_aligned, 100*Pseudo_Prop_masked, linewidth = 4, color = color_list[c_ind], label = "%s"%countries_labels[c_ind])
                    
                #ax_twin.scatter(inds_dates, 100*Pseudo_Prop_masked, marker = ".", color = color_list[k])
    
                ### Separate figure for relative fitness vs change in proportion
                gamma_prop_masked = ma.masked_array(gamma_prop, mask = SI_mask)
                
                if c_ind == 0:
                    ### only plot the fitness of the main country investigated
                    ax2.set_title(lineage_list[k].replace("/", " + ") + " %s"%countries_labels[0][:3].upper(), fontsize = 32)
                    ax2_twin.plot(inds_dates, gamma_prop_masked, color = color_list[c_ind], linewidth = 4, label = "%s"%countries_labels[c_ind])
                    ax2.fill_between(inds_dates_gamma, gamma_SI_min, gamma_SI_max, color = color_list[c_ind], alpha = 0.3, label = "%s"%countries_labels[c_ind])
                else:
                    ### placeholder for legend
                    ax2.fill_between(inds_dates_gamma, gamma_SI_min, gamma_SI_max, color = color_list[c_ind], alpha = 0.3, label = "%s"%countries_labels[c_ind])
                    ax2_twin.plot(inds_dates, gamma_prop_masked, color = color_list[c_ind], linewidth = 4, label = "%s"%countries_labels[c_ind])
                
                #ax2_twin.scatter(inds_dates, gamma_prop_masked, marker = ".", color = "orange")
                            
                status_list.append(lab_status)
                status_list_pseudo.append(lab_status_pseudo)
                save_plot.append(True)
            else:
                print("No lineage in group %s has E[Susceptible] available, if needed, first compute it in main config"%lineage_list[k], ES_df_dir)
                status_list.append("%s : No data"%countries_labels[c_ind])
                status_list_pseudo.append("%s : No data"%countries_labels[c_ind])
                save_plot.append(False)
            
            lineage_list_countries.append("%s %s"%(lineage_list[k].replace("/", "+"), countries_labels[0][:3].upper()))
               
             
        if np.any(save_plot):
            if len(other_props) != 0:
                
                sums = [np.sum(other_props[ci][:60]) for ci in range(len(other_props))]
                test = np.array(sums)>np.sum(Prop_ref[:60])
                indx_list = np.where(test)[0]
                    
                if len(indx_list) > 0:
                    indx = np.argmax(sums)
                else:
                    indx = None
                    
                if indx is not None:
                    inds_dates_gamma = other_t[indx]
                    gamma_SI_min, gamma_SI_max = other_trends[indx]
                    ax_twin.fill_between(inds_dates_gamma, gamma_SI_min, gamma_SI_max, color = other_cols[indx], alpha = 0.3, label = "%s"%other_labs[indx])
                    ax.plot(other_t_all[indx], 100*other_prop_all[indx], linewidth = 6, color = other_cols[indx], label = "%s"%other_labs[indx])
                    lim_0, lim_1 = ax_twin.get_ylim()
                    
                sorted_inds = np.argsort(sums)[::-1]
                lw = np.linspace(3, 5, len(other_cols))[::-1]
                alph = np.linspace(0.3, 0.6, len(other_cols))
                
                for cs in range(len(sorted_inds)):
                    ci = sorted_inds[cs]
                    plot = True
                    if indx is not None:
                        if ci == indx:
                            plot = False     
                    if plot:
                        gamma_SI_min, gamma_SI_max = other_trends[ci]
                        ax_twin.fill_between(other_t[ci], gamma_SI_min, gamma_SI_max, color = other_cols[ci], alpha = 0., label = "%s"%other_labs[ci])
                        ax.plot(other_t_all[ci], 100*other_prop_all[ci], linewidth = lw[cs], alpha = alph[cs] ,color = other_cols[ci], label = "%s"%other_labs[ci])

                if indx is not None:
                    ax.plot(other_t_all[indx], 100*other_prop_all[indx], linewidth = 6, color = other_cols[indx])
                
                ax.plot(t_prop_aligned_0, 100*Pseudo_Prop_masked_0, linewidth = 6, color = color_list[0])                
            
            ax_twin.axhline(xmin = 0, xmax = len(all_dates), ls = "--", linewidth = 2, color = "black")
            
            #ymin1, ymax1 = ax.get_ylim()
            ymin2, ymax2 = lim_0, lim_1
            #ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            #ax.set_ylim((ymin1, ymax1)
            ax_twin.set_ylim((ymin2, ymax2))
            #loc0 = min(np.abs(ymin1)/(np.abs(ymin1)+np.abs(ymax1)), np.abs(ymax1)/(np.abs(ymin1)+np.abs(ymax1)))
            #mpl_axes_aligner.align.yaxes(ax, 0, ax_twin, 0, loc0)
        
            
            ax.set_xticks(perday)
            ax.set_xticklabels(date_ticks,
                rotation = 45, horizontalalignment = "right")
            
            if (x_min is not None):
                ax.set_xlim((x_min, x_max))
                ax_twin.set_xlim((x_min, x_max))
            
            
            lab_k_fn = (lab_k.replace("/", "_")).replace("*","").replace("+", "_") ## for filename 
            if len(lab_k_fn) > 10: # can't be too long
                lab_k_fn = lab_k_fn[:10] + "_et_al"
            
            ax_twin.set_zorder(-1) ### send the legend of ax_twin in the background
            ax_twin.patch.set_visible(False) ### make sure that the axis plots do not get hiden in the background
            ax.patch.set_visible(False) ### make sure that the axis plots ddo not get hiden in the background
            
            ax_twin.set_ylabel("Relative fitness", fontsize = 20)
            ax.set_ylabel("Lineage Frequency (daily %)", fontsize = 20)
            ax.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = np.ceil(len(lineage_list)/12).astype(int))
            ax_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = np.ceil(len(lineage_list)/12).astype(int))
            pdf = PdfPages(sys.argv[w_save]+"/relative_fitness_prop_%s.pdf"%lineage_list[k].replace("/", "_"))
            pdf.savefig(fig, bbox_inches = "tight")
            pdf.close()
            fig.savefig(sys.argv[w_save]+"/relative_fitness_prop_%s.svg"%lineage_list[k].replace("/", "_"), bbox_inches = "tight")
            plt.close()
            
            ### vs proportions
            ax2.set_xticks(perday)
            ax2.set_xticklabels(date_ticks,
                rotation = 45, horizontalalignment = "right")
            
            
            if x_min is not None:
                ax2.set_xlim((x_min, min(x_max, x_max)))
                ax2_twin.set_xlim((x_min, min(x_max, x_max)))
            
            ymin1, ymax1 = ax2.get_ylim()
            ymin2, ymax2 = ax2_twin.get_ylim()
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            
            #loc0 = min(np.abs(ymin1)/(np.abs(ymin1)+np.abs(ymax1)), np.abs(ymax1)/(np.abs(ymin1)+np.abs(ymax1)))
            #mpl_axes_aligner.align.yaxes(ax2, 0, ax2_twin, 0, loc0)
            mpl_axes_aligner.align.yaxes(ax2, 0, ax2_twin, 0, 0.5)
            
            if (ymin1/ymin2 >0.5) or (ymax1/ymax2>0.5) or (ymin2/ymin1 >0.5) or (ymax2/ymax1>0.5):
                ax2.set_ylim((ymin, ymax))
                ax2_twin.set_ylim((ymin, ymax))   
            
                
            ax2.axhline(xmin = 0, xmax = len(day_prop), ls = "--", linewidth = 2, color = "black")
            ax2.set_ylabel("Relative fitness $\gamma_y$", fontsize = 20)
            ax2_twin.set_ylabel("Change in proportion $\gamma_{prop}$", fontsize = 20)
            ax2.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = np.ceil(len(lineage_list)/12).astype(int))
            ax2_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = np.ceil(len(lineage_list)/12).astype(int))
            pdf_2 = PdfPages(sys.argv[w_save]+"/relative_fitness_%s.pdf"%lineage_list[k].replace("/", "_"))
            pdf_2.savefig(fig2, bbox_inches = "tight")
            pdf_2.close()
            fig2.savefig(sys.argv[w_save]+"/relative_fitness_%s.svg"%lineage_list[k].replace("/", "_"), bbox_inches = "tight")
            plt.close()
        
    return status_list, status_list_pseudo, lineage_list_countries

num_groups = int(sys.argv[7])
w_save = 8
k = 9
lineage_list = []
color_list = []
custom_col = sns.color_palette("Set2", 100) 
s = 0

countries_dir_list = str(sys.argv[1]).split("/new/")
countries_list = str(sys.argv[2]).split("/")
countries_labels = str(sys.argv[3]).split("/")

for i in range(num_groups):
    lineage_list.append(str(sys.argv[k+i]))
    s += 1

if "/new/" in str(sys.argv[s+k]):
    cols = str(sys.argv[s+k]).split("/new/")
    for c in cols:
        try:
            if "/" not in c:
                color_list.append(c)
            else:
                split_col = str(c).split("/")
                color_list.append(tuple([float(split_col[c]) for c in range(len(split_col))])) ### anything else is error
        except:
            rand_num = np.random.choice(1, 100)
            if s<len(custom_col):
                color_list.append(custom_col[s])
            else:
                color_list.append(sns.color_palette("rocked", rand_num)[0])
            s +=1
else:
    for j in range(len(countries_list)):
        try:
            if "/" not in str(sys.argv[s+k+j]):
                color_list.append(str(sys.argv[s+k+j]))
            else:
                split_col = str(sys.argv[s+k+j]).split("/")
                color_list.append(tuple([float(split_col[c]) for c in range(len(split_col))])) ### anything else is error
        except:
            rand_num = np.random.choice(1, 100)
            if s<len(custom_col):
                color_list.append(custom_col[s])
            else:
                color_list.append(sns.color_palette("rocked", rand_num)[0])
            s +=1

status_list, status_list_pseudo, lineage_list_countries = plot_fit(countries_dir_list, countries_list, countries_labels, lineage_list, color_list, w_save)
status = pd.DataFrame({"lineage":lineage_list_countries, "Relevant lineages":status_list, "Pseudo Grouping": status_list_pseudo})
status.to_csv(sys.argv[w_save]+"/plot_status.csv")
 

        