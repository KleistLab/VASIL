#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:53:58 2023

@author: raharinirina
"""
import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import pdb
import seaborn as sns
import pickle

def moving_average(X, window = 7):
    u = np.zeros(len(X))
    u[:window] = X[:window]
    for i in range(window, len(X)):
        u[i] = np.mean(X[i-window:i+1])
    
    return u

ES_lin_dir = sys.argv[1]
S_mean_file = sys.argv[2]
lineage_freq = pd.read_csv(sys.argv[3])
threshold = float(sys.argv[4])

#ES_df = pd.read_csv("demo/results/Immunological_Landscape/Susceptible_SpikeGroup_lineage_XXX_all_PK.csv")
#lineage_freq = pd.read_csv("demo/results/Daily_Lineages_Freq.csv")
#threshold = 0 #(percent)
#variant = "BA.5.1"

# needs to be updated to allow individual weighting 
S_mean_df = pd.read_csv(S_mean_file)
S_all_mean = S_mean_df.to_numpy()[:, S_mean_df.columns != "Days"].astype(float)
t_dates = S_mean_df["Days"]

# processing of frequency data
try:
    lineage_freq.drop(columns = "Unnamed: 0", inplace = True)
except:
    pass

lineage_freq = lineage_freq#[lineage_freq['date'].isin(t_dates)]
freqs = lineage_freq.loc[:, lineage_freq.columns != 'date']

# imputing frequencies below threshold and normalization
freqs = freqs.mask(freqs < threshold)
freqs = freqs.fillna(0)
col_sums = freqs.sum(axis = 1).values
freqs = freqs.divide(col_sums, axis="rows")
freqs = freqs.fillna(0)
lineage_freq.loc[:, lineage_freq.columns != 'date'] = freqs
day_prop = lineage_freq["date"].tolist()
t_prop = np.arange(len(day_prop)).astype(int)

import matplotlib
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)

file = open("Spikegroups_membership.pck", "rb")
Pseudogroup_dic = pickle.load(file)
file.close()    
def plot_fit(ES_df_dir, lineage_list, color_list, w_save = len(sys.argv)-1):
    # plotting
    PreFig(xsize = 20, ysize = 20)
    fig = plt.figure(figsize = (9, 7))
    ax = fig.add_subplot(1, 1, 1)
    ### end of observation line
    ax.axvline(x = len(t_dates) - 1, ymin = -1, ymax = 1, ls = "--", linewidth = 2, color = "grey")

    # different axis for proportions
    ax_twin = ax.twinx()

    status_list = []
    
    Pseudo_prop = np.zeros((len(t_prop)))
    for k in range(len(lineage_list)):
        splited_var = np.array(lineage_list[k].split("/"))
        splited_var = splited_var[~(splited_var == "")]
        splited_var = splited_var[~(splited_var == " ")]
        num_avail = 0
        
        # plotting
        PreFig(xsize = 20, ysize = 20)
        fig_k = plt.figure(figsize = (9, 7))
        ax_k = fig_k.add_subplot(1, 1, 1)
        ### end of observation line
        ax_k.axvline(x = len(t_dates) - 1, ymin = -1, ymax = 1, ls = "--", linewidth = 2, color = "grey")

        # different axis for proportions
        ax_k_twin = ax_k.twinx()
        
        for x in range(len(splited_var)):
            lineage = splited_var[x]
            try:
                ES_df = pd.read_csv(ES_df_dir+"/Susceptible_SpikeGroup_%s_all_PK.csv"%lineage)
                num_avail +=1
            except:
                try:
                    ES_df = pd.read_csv("results/Immunological_Landscape_ALL/Susceptible_SpikeGroup_%s_all_PK.csv"%lineage)
                    num_avail +=1
                except:
                    try:
                        ES_df = pd.read_csv("results/Immunological_Landscape/Susceptible_SpikeGroup_%s_all_PK.csv"%lineage)
                        num_avail +=1
                    except:
                        print("Computation needed: Excpected Susceptible file is not available for %s"%lineage)
            # processing of susceptibles 
            try:
                ES_df.drop(columns = "Unnamed: 0", inplace = True)
            except:
                pass
            
            es_cols = ES_df.columns
            ES_df = ES_df[ES_df['Days'].isin(t_dates)]
            if x == 0:
                ES_sum = ES_df.to_numpy()[:, es_cols!="Days"].astype(float)
            else:
                ES_sum += ES_df.to_numpy()[:, es_cols!="Days"].astype(float)
            
            # change in relative frequency from genomic surveillance data 
            Pseudo_prop = np.zeros((len(t_prop)))
            if lineage in list(Pseudogroup_dic.keys()):
                if "Spike. " + Pseudogroup_dic[lineage] in lineage_freq.columns.astype(str):
                    if x == 0:
                        Pseudo_Prop = moving_average(lineage_freq["Spike. " + Pseudogroup_dic[lineage]], window = 14)
                    else:
                        Pseudo_Prop += moving_average(lineage_freq["Spike. " + Pseudogroup_dic[lineage]], window = 14)
    
                    #Pseudo_Prop[Pseudo_Prop < threshold] = 0        
                    #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
                elif lineage in lineage_freq.columns.astype(str):
                    if x == 0:
                        Pseudo_Prop = moving_average(lineage_freq[lineage], window = 14)
                    else:
                        Pseudo_Prop += moving_average(lineage_freq[lineage], window = 14)
                    #Pseudo_Prop[Pseudo_Prop < threshold] = 0
                    #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
            else:
                if x == 0:
                    Pseudo_Prop = ma.masked_array(np.zeros(len(t_prop)), mask = np.ones(len(t_prop), dtype = bool))
                else:
                    Pseudo_Prop += ma.masked_array(np.zeros(len(t_prop)), mask = np.ones(len(t_prop), dtype = bool))
        
        
        if num_avail !=0:
            ES_ranges = ES_sum/num_avail# compute the mean
            Pseudo_prop = Pseudo_prop/num_avail 
        else:
            print("Error: There are no E[Susceptible] files for any lineage in %s"%lineage)
            
        """
        # calculation of change in relative frequency from model
        gamma_prop = np.zeros(len(t_dates))
        for l in range(len(t_dates)-1):
            if Pseudo_Prop[l] == 0 or Pseudo_Prop[l+1] == 0:
                gamma_prop[l] = float('nan')
            else:
                gamma_prop[l] = Pseudo_Prop[l+1]/Pseudo_Prop[l] -1
        ax_twin.plot(t_dates, gamma_prop, color = color_list[k], label = lineage)
        """   
        
        gamma_SI = np.zeros((len(t_dates), ES_ranges.shape[1]))
        
        for i in range(ES_ranges.shape[1]):
            S_x = ES_ranges[:, i]
            S_mean = S_all_mean[:, i]
            
            gamma_SI[:, i] = np.divide(S_x - S_mean, S_mean, out = S_x, where = S_mean != 0)
        
        # get min max gamma over PK at each timepoints
        gamma_SI_min, gamma_SI_max = np.min(gamma_SI, axis = 1), np.max(gamma_SI, axis = 1)
        
        inds_dates = np.arange(0,len(t_dates),1)
        ax.fill_between(inds_dates, gamma_SI_min, gamma_SI_max, color = color_list[k], alpha = 0.3, label = lineage_list[k])
        ax_twin.plot(t_prop, Pseudo_Prop, linewidth = 3, color = color_list[k], label = lineage_list[k])
        status_list.append("Done")
        
        ax_k.fill_between(inds_dates, gamma_SI_min, gamma_SI_max, color = color_list[k], alpha = 0.3, label = lineage_list[k])
        ax_k_twin.plot(t_prop, Pseudo_Prop, linewidth = 3, color = color_list[k], label = lineage_list[k])
        ax_k.axhline(xmin = 0, xmax = len(t_dates), ls = "--", linewidth = 2, color = "black")
        
        ymin1, ymax1 = ax_k.get_ylim()
        ymin2, ymax2 = ax_k_twin.get_ylim()
        ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
        ax_k.set_ylim((ymin, ymax))
        ax_k_twin.set_ylim((ymin, ymax))
        try:
            x_min = list(t_dates).index(str(sys.argv[5]))
            #x_max = list(t_dates).index(str(sys.argv[6]))
            x_min1 = day_prop.index(str(sys.argv[5]))
            x_max1 = day_prop.index(str(sys.argv[6]))
            x_max = x_max1
        except:
            x_min = None
    
        if (x_min is not None):
            ax_k.set_xlim((x_min, x_max))
            ax_k_twin.set_xlim((x_min1, x_max1))
            t_dates_show = np.array(t_dates)[x_min:x_max+1]
        else:
            t_dates_show = t_dates
    
        if len(t_dates_show)>200:
            pp = 7*4
        else:
            pp = min(len(t_dates_show), 14)
        
        perday = np.arange(0,len(t_dates_show), pp)
        date_ticks = t_dates_show[perday].tolist()
        if t_dates[len(t_dates) - 1] not in date_ticks:
            n=list(t_dates).index(date_ticks[-1])+pp
            while n<len(t_dates)-1:
                date_ticks.append(t_dates[n])
                perday = np.append(perday, n)
                n += pp
            date_ticks.append(t_dates[len(t_dates) - 1])
            perday = np.append(perday, len(t_dates) - 1)
        
        change = len(date_ticks)
        
        if day_prop[len(day_prop) - 1] not in date_ticks:
            n=list(day_prop).index(date_ticks[-1])+pp
            while n<len(day_prop)-1:
                date_ticks.append(day_prop[n])
                perday = np.append(perday, n)
                n += pp
            date_ticks.append(day_prop[len(t_prop) - 1])
            perday = np.append(perday, len(t_prop) - 1)
    
        if x_min is not None:
            perday_orig = []
            for i in range(len(np.array(date_ticks)[:change])):
                perday_orig.append(list(t_dates).index(date_ticks[i]))
            try:
                for j in range(len(np.array(date_ticks[change:]))):
                    perday_orig.append(list(day_prop).index(date_ticks[change+j]))
            except:
                pass
        else:
            perday_orig = perday
            
        ax_k.set_xticks(perday_orig)
        ax_k.set_xticklabels(date_ticks,
            rotation = 45, horizontalalignment = "right")
    
        #ax_twin.set_ylim((-0.02, 0.02))
        
        ax_k.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
        ax_k_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
        ax_k.set_ylabel("Relative fitness", fontsize = 20)
        ax_k_twin.set_ylabel("Variant abundance (daily)", fontsize = 20)
        pdf_k = PdfPages(sys.argv[w_save]+"/relative_fitness_%s.pdf"%lineage_list[k].replace("/", "_"))
        pdf_k.savefig(fig_k, bbox_inches = "tight")
        pdf_k.close()
 
        fig_k.savefig(sys.argv[w_save]+"/relative_fitness_%s.svg"%(lineage_list[k].replace("/", "_")), bbox_inches = "tight")
        plt.close()
        
    ax.axhline(xmin = 0, xmax = len(t_dates), ls = "--", linewidth = 2, color = "black")
    
    ymin1, ymax1 = ax.get_ylim()
    ymin2, ymax2 = ax_twin.get_ylim()
    ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
    ax.set_ylim((ymin, ymax))
    ax_twin.set_ylim((ymin, ymax))
    try:
        x_min = list(t_dates).index(str(sys.argv[5]))
        #x_max = list(t_dates).index(str(sys.argv[6]))
        x_min1 = day_prop.index(str(sys.argv[5]))
        x_max1 = day_prop.index(str(sys.argv[6]))
        x_max = x_max1
    except:
        x_min = None

    if (x_min is not None):
        ax.set_xlim((x_min, x_max))
        ax_twin.set_xlim((x_min1, x_max1))
        t_dates_show = np.array(t_dates)[x_min:x_max+1]
    else:
        t_dates_show = t_dates

    if len(t_dates_show)>200:
        pp = 7*4
    else:
        pp = min(len(t_dates_show), 14)
    
    perday = np.arange(0,len(t_dates_show), pp)
    date_ticks = t_dates_show[perday].tolist()
    if t_dates[len(t_dates) - 1] not in date_ticks:
        n=list(t_dates).index(date_ticks[-1])+pp
        while n<len(t_dates)-1:
            date_ticks.append(t_dates[n])
            perday = np.append(perday, n)
            n += pp
        date_ticks.append(t_dates[len(t_dates) - 1])
        perday = np.append(perday, len(t_dates) - 1)
    
    change = len(date_ticks)
    
    if day_prop[len(day_prop) - 1] not in date_ticks:
        n=list(day_prop).index(date_ticks[-1])+pp
        while n<len(day_prop)-1:
            date_ticks.append(day_prop[n])
            perday = np.append(perday, n)
            n += pp
        date_ticks.append(day_prop[len(t_prop) - 1])
        perday = np.append(perday, len(t_prop) - 1)

    if x_min is not None:
        perday_orig = []
        for i in range(len(np.array(date_ticks)[:change])):
            perday_orig.append(list(t_dates).index(date_ticks[i]))
        try:
            for j in range(len(np.array(date_ticks[change:]))):
                perday_orig.append(list(day_prop).index(date_ticks[change+j]))
        except:
            pass
    else:
        perday_orig = perday
        
    ax.set_xticks(perday_orig)
    ax.set_xticklabels(date_ticks,
        rotation = 45, horizontalalignment = "right")

    #ax_twin.set_ylim((-0.02, 0.02))
    
    ax.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
    ax_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
    ax.set_ylabel("Relative fitness", fontsize = 20)
    ax_twin.set_ylabel("Variant abundance (daily)", fontsize = 20)
    pdf = PdfPages(sys.argv[w_save]+"/relative_fitness_groups.pdf")
    pdf.savefig(fig, bbox_inches = "tight")
    pdf.close()
 
    fig.savefig(sys.argv[w_save]+"/relative_fitness_groups.svg", bbox_inches = "tight")
    plt.close()
    return status_list

num_groups = int(sys.argv[7])
w_save = 8
k = 9
lineage_list = []
color_list = []
custom_col = sns.color_palette("Set2", 100) 
s = 0
for i in range(num_groups):
    lineage_list.append(str(sys.argv[k+i]))
    try:
        color_list.append(str(sys.argv[k+num_groups+i]))
    except:
        rand_num = np.random.choice(1, 100)
        if s<len(custom_col):
            color_list.append(custom_col[s])
        else:
            color_list.append(sns.color_palette("rocked", rand_num)[0])
        s +=1

status_list = plot_fit(ES_lin_dir, lineage_list, color_list, w_save)
    
status = pd.DataFrame({"lineage":lineage_list, "relative_advantage":status_list})
status.to_csv(sys.argv[w_save]+"/plot_status.csv")
 

        