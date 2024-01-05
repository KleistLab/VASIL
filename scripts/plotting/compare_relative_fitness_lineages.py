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
Requires system arguments to compare several trends from different simulations

sys.argv[1]: "path/to/Trend_1/new/path/to/Trend_2/new/path/to/Trend_3/..." ### as many Trends as desired but always separated by /new/ and folder must include all Trend results
sys.argv[2]: "Trend_1/Trend_2/Trend_3/..." ### Trend names separated by / 
sys.argv[2]: "Trend_1/Trend_2/Trend_3/..." ### Trend labels separated by /
sys.argv[4]: percentage to filter spikegroups #e.g. 0.01 
sys.argv[5]: dates to start the comparison for each variants, separated by /, if only one then they will be all the same ### e.g. "2022-04-09/2022-04-26" for 2 variants
sys.argv[6]: dates to end the comparison for each variants, separated by /, if only one then they will be all the same ### e.g. "2022-08-21/2023-04-15" for 2 variants
sys.argv[7]: integer number of variants to check
sys.argv[8]: folder path to store the results 
sys.argv[9:9+sys.argv[7]]: list all lineages to simulate ### e.g "BE.10.ALL", "BE.10.ALL" 
sys.argv[9+sys.argv[7]:len(sys.argv)]: list 1 color for each Trends by their order in sys.argv[2] ### e.g. "red" "orange" for two Trends
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

   
def plot_fit(Trends_dir_list, Trends_subdir_list, Trends_labels, lineage_list, color_list, w_save = len(sys.argv)-1):
    status_list = []
    lineage_list_Trends = []    
    for k in range(len(lineage_list)):
        t_dates_list = []
        S_all_mean_list = []
        all_dates = []
        all_prop_dates = []
        lineage_freqs_list = []
        prop_mask_list = []
        day_prop_list = []
        for c_ind in range(len(Trends_subdir_list)):
            if "ImL_ALL_vs_Vacc_ver2" in Trends_subdir_list[c_ind]:
                S_mean_file = Trends_dir_list[c_ind]+"/Susceptible_weighted_mean_over_spikegroups_vs_Vacc_ver2_all_PK.csv"
            elif "ImL_ALL_vs_Vacc_ver1" in Trends_subdir_list[c_ind]:
                S_mean_file = Trends_dir_list[c_ind]+"/Susceptible_weighted_mean_over_spikegroups_vs_Vacc_ver1_all_PK.csv"
            else:
                S_mean_file = Trends_dir_list[c_ind]+"/Susceptible_weighted_mean_over_spikegroups_all_PK.csv"
            
            print("-------------------------------------------------------------")
            print(Trends_labels[c_ind], "open S_mean", S_mean_file)
            # needs to be updated to allow individual weighting 
            S_mean_df = pd.read_csv(S_mean_file)
            S_all_mean_list.append(S_mean_df.to_numpy()[:, S_mean_df.columns != "Days"].astype(float))
            t_dates_list.append(S_mean_df["Days"].tolist())
            all_dates += list(S_mean_df["Days"])
            
            # processing of frequency data
            file = open(Trends_dir_list[c_ind].replace("/vaccination", "")+"/Spikegroups_membership.pck", "rb") #hard-coded for vaccination simulations
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
        
        t_prop_all = np.arange(len(all_prop_dates)).astype(int)
        inds_dates_all = np.arange(len(all_dates)).astype(int)
        
        ax.axvline(x = len(all_dates) - 1, ymin = -1, ymax = 1, ls = "--", linewidth = 2, color = "grey")
        # different axis for proportions
        ax_twin = ax.twinx()
        
        ### Separate figure for relative fitness vs change in proportion
        PreFig(xsize = 20, ysize = 20)
        fig2 = plt.figure(figsize = (15, 7))
        ax2 = fig2.add_subplot(1, 1, 1)
        # different axis for proportions
        ax2_twin = ax2.twinx()
        
        for c_ind in range(len(Trends_subdir_list)):
            file = open(Trends_dir_list[c_ind].replace("/vaccination", "")+"/Spikegroups_membership.pck", "rb") #hard-coded for vaccination simulations
            Pseudogroup_dic = pickle.load(file)
            file.close() 
            
            ES_df_dir = Trends_dir_list[c_ind]+"/"+Trends_subdir_list[c_ind]
            day_prop = day_prop_list[c_ind]
            t_prop = t_prop_all[list(all_prop_dates).index(day_prop[0]):list(all_prop_dates).index(day_prop[-1])+1]
            prop_mask = prop_mask_list[c_ind]
            t_dates = t_dates_list[c_ind]
            inds_dates = inds_dates_all[list(all_dates).index(t_dates[0]):list(all_dates).index(t_dates[-1])+1]
            S_all_mean = S_all_mean_list[c_ind]
            variant_freq = lineage_freqs_list[c_ind].copy()
            Pseudo_done_global = []
        
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
                        print("-------------------------------------------------------------")
                        print(Trends_labels[c_ind], ES_df_dir+"/Susceptible_SpikeGroup_%s_all_PK.csv"%Pseudogroup_dic[lineage])
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
                                    print("Computation needed: Expected Susceptible file is not available for %s"%lineage)
                                    run = False
                                    
                # processing of Proportions data
                if run:   
                   # change in relative frequency from genomic surveillance data 
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
            
            lab_status = "%s : "%Trends_labels[c_ind] + lab_status
                
                
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
                prop_mask_aligned = np.zeros(len(t_dates))
                gamma_prop = np.zeros(len(t_dates))
                SI_mask = np.zeros(len(t_dates)).astype(bool)
                for l in range(len(t_dates)):
                    #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
                    ## Dates are already set to always be to be successive
                    if t_dates[l] in day_prop:
                        w_l = list(day_prop).index(t_dates[l])
                        Pseudo_Prop_aligned[l] = Pseudo_Prop[w_l]
                        prop_mask_aligned[l] = prop_mask.tolist()[w_l]
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
                
                # calculation of relative fitness
                gamma_SI = np.zeros((len(t_dates), ES_ranges.shape[1]))
                for i in range(ES_ranges.shape[1]):
                    S_x = ES_ranges[:, i]
                    S_mean = S_all_mean[:, i]
                    gamma_SI[:, i] = np.divide(S_x - S_mean, S_mean, out = S_x, where = (S_mean != 0)&(S_x!=0))
                
                # get min max gamma over PK at each timepoints
                SI_mask = np.array(SI_mask) + prop_mask[:len(inds_dates)] ### props are already aligned with indicence date
                gamma_SI_min, gamma_SI_max = np.min(gamma_SI, axis = 1), np.max(gamma_SI, axis = 1)
                gamma_SI_max = ma.masked_array(gamma_SI_max, mask = SI_mask)
                gamma_SI_min = ma.masked_array(gamma_SI_min, mask = SI_mask)
                Pseudo_Prop_masked = ma.masked_array(Pseudo_Prop_aligned, mask = prop_mask_aligned)
    
                ax_twin.fill_between(inds_dates, gamma_SI_min, gamma_SI_max, color = color_list[c_ind], alpha = 0.3, label = lab_k + " -- %s"%Trends_labels[c_ind])
                
                if ("ImL_ALL_vs_Vacc_ver2" not in Trends_subdir_list[c_ind]) and ("ImL_ALL_vs_Vacc_ver1" not in Trends_subdir_list[c_ind]):
                    ax.plot(inds_dates, 100*Pseudo_Prop_masked, linewidth = 4, color = color_list[c_ind], label = lab_k + " -- %s"%Trends_labels[c_ind])
                else:
                    ax.plot(inds_dates, 100*Pseudo_Prop_masked, linewidth = 4, color = color_list[c_ind],  alpha = 0., label = lab_k + " -- %s"%Trends_labels[c_ind])    
                #ax_twin.scatter(t_prop, Pseudo_Prop_masked, marker = ".", color = color_list[k])
    
                ### Separate figure for relative fitness vs change in proportion
                gamma_prop_masked = ma.masked_array(gamma_prop, mask = SI_mask)
                ax2.fill_between(inds_dates, gamma_SI_min, gamma_SI_max, color = color_list[c_ind], alpha = 0.3, label = lab_k + " -- %s"%Trends_labels[c_ind])
                
                if ("ImL_ALL_vs_Vacc_ver2" not in Trends_subdir_list[c_ind]) and ("ImL_ALL_vs_Vacc_ver1" not in Trends_subdir_list[c_ind]):
                    ax2_twin.plot(inds_dates, gamma_prop_masked, color = color_list[c_ind], linewidth = 4,  label=lab_k + " -- %s"%Trends_labels[c_ind])
                else:
                    ax2_twin.plot(inds_dates, gamma_prop_masked, color = color_list[c_ind], linewidth = 4,  alpha = 0., label=lab_k + " -- %s"%Trends_labels[c_ind])
                #ax2_twin.scatter(inds_dates, gamma_prop_masked, marker = ".", color = "orange")
                            
                status_list.append(lab_status)
            
            else:
                print("No lineage in group %s has E[Susceptible] available, if needed, first compute it in main config"%lineage_list[k])
                status_list.append("%s : No data"%Trends_labels[c_ind])
            
            lineage_list_Trends.append("%s"%(lineage_list[k]))
        
        ax_twin.axhline(xmin = 0, xmax = len(all_dates), ls = "--", linewidth = 2, color = "black")
        
        ymin1, ymax1 = ax.get_ylim()
        ymin2, ymax2 = ax_twin.get_ylim()
        #ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
        ax.set_ylim((ymin1, ymax1))
        #ax_twin.set_ylim((ymin2, ymax2))
        #loc0 = min(np.abs(ymin1)/(np.abs(ymin1)+np.abs(ymax1)), np.abs(ymax1)/(np.abs(ymin1)+np.abs(ymax1)))
        #mpl_axes_aligner.align.yaxes(ax, 0, ax_twin, 0, loc0)
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
                x_min = 0
                
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

        ax_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = 1)
        ax.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = 1)
        ax_twin.set_ylabel("Relative fitness", fontsize = 20)
        ax.set_ylabel("Spikegroup Frequency (daily %)", fontsize = 20)
        pdf = PdfPages(sys.argv[w_save]+"/relative_fitness_prop_%s.pdf"%lab_k_fn)
        pdf.savefig(fig, bbox_inches = "tight")
        pdf.close()
        fig.savefig(sys.argv[w_save]+"/relative_fitness_prop_%s.svg"%lab_k_fn, bbox_inches = "tight")
        plt.close()
        
        ### vs proportions
        ax2.set_xticks(perday)
        ax2.set_xticklabels(date_ticks,
            rotation = 45, horizontalalignment = "right")
        
        if x_min is not None:
            ax2.set_xlim((x_min, min(x_max, inds_dates[-1])))
            ax2_twin.set_xlim((x_min, min(x_max, inds_dates[-1])))
        
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
        ax2.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = 1)
        ax2_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = 1)
        pdf_2 = PdfPages(sys.argv[w_save]+"/relative_fitness_%s.pdf"%lab_k_fn)
        pdf_2.savefig(fig2, bbox_inches = "tight")
        pdf_2.close()
        fig2.savefig(sys.argv[w_save]+"/relative_fitness_%s.svg"%lab_k_fn, bbox_inches = "tight")
        plt.close()
        
    return status_list, lineage_list_Trends

num_groups = int(sys.argv[7])
w_save = 8
k = 9
lineage_list = []
color_list = []
custom_col = sns.color_palette("Set2", 100) 
s = 0

Trends_dir_list = str(sys.argv[1]).split("/new/")
Trends_subdir_list = str(sys.argv[2]).split("/new/")
Trends_labels = str(sys.argv[3]).split("/")

for i in range(num_groups):
    lineage_list.append(str(sys.argv[k+i]))
    s+=1


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
    for j in range(len(Trends_subdir_list)):
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


status_list, lineage_list_Trends = plot_fit(Trends_dir_list, Trends_subdir_list, Trends_labels, lineage_list, color_list, w_save)

status = pd.DataFrame({"lineage":lineage_list_Trends, "spikegroups_found":status_list})
status.to_csv(sys.argv[w_save]+"/plot_status.csv")
 

        