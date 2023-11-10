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

lineage_freq = lineage_freq[lineage_freq['date'].isin(t_dates)]
freqs = lineage_freq.loc[:, lineage_freq.columns != 'date']

# imputing frequencies below threshold and normalization
freqs = freqs.mask(freqs < threshold)
freqs = freqs.fillna(0)
col_sums = freqs.sum(axis = 1).values
freqs = freqs.divide(col_sums, axis="rows")
freqs = freqs.fillna(0)
lineage_freq.loc[:, lineage_freq.columns != 'date'] = freqs

import matplotlib
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)
    
def plot_fit(ES_df_list, lineage_list, color_list, w_save = len(sys.argv)-1):
    # plotting
    PreFig(xsize = 20, ysize = 20)
    fig = plt.figure(figsize = (15, 7))
    ax = fig.add_subplot(1, 1, 1)
    ### end of observation line
    ax.axvline(x = len(t_dates) - 1, ymin = -1, ymax = 1, ls = "--", linewidth = 2, color = "grey")
    #ax.set_ylim((-0.02, 0.02))

    # different axis for proportions
    ax_twin = ax.twinx()
    status_list = []
    for k in range(len(ES_df_list)):
        ES_df = ES_df_list[k]
        lineage = lineage_list[k]
        # processing of susceptibles 
        try:
            ES_df.drop(columns = "Unnamed: 0", inplace = True)
        except:
            pass
        
        es_cols = ES_df.columns
        ES_df = ES_df[ES_df['Days'].isin(t_dates)]
        ES_ranges = ES_df.to_numpy()[:, es_cols!="Days"].astype(float)
        
        # calculation of change in relative frequency from model
        gamma_SI = np.zeros((len(t_dates), ES_ranges.shape[1]))
        
        for i in range(ES_ranges.shape[1]):
            S_x = ES_ranges[:, i]
            S_mean = S_all_mean[:, i]
        
            gamma_SI[:, i] = np.divide(S_x - S_mean, S_mean, out = S_x, where = S_mean != 0)
        
        # get min max gamma over PK at each timepoints
        gamma_SI_min, gamma_SI_max = np.min(gamma_SI, axis = 1), np.max(gamma_SI, axis = 1)
        
        # change in relative frequency from genomic surveillance data 
        if "Spike. " + lineage in lineage_freq.columns.astype(str):
            Pseudo_Prop = moving_average(lineage_freq["Spike. " + lineage], window = 14)
            #Pseudo_Prop[Pseudo_Prop < threshold] = 0        
            #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
        elif lineage in lineage_freq.columns.astype(str):
            Pseudo_Prop = moving_average(lineage_freq[lineage], window = 14)
            #Pseudo_Prop[Pseudo_Prop < threshold] = 0
            #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
        else:
            Pseudo_Prop = np.zeros(len(t_dates))
        
        gamma_prop = np.zeros(len(t_dates))
        for l in range(len(t_dates)-1):
            if Pseudo_Prop[l] == 0 or Pseudo_Prop[l+1] == 0:
                gamma_prop[l] = float('nan')
            else:
                gamma_prop[l] = Pseudo_Prop[l+1]/Pseudo_Prop[l] -1
    
        ax.fill_between(t_dates, gamma_SI_min, gamma_SI_max, color = color_list[k], alpha = 0.3, label = lineage)
        ax_twin.plot(t_dates, gamma_prop, color = color_list[k])
        status_list.append("Done")

    #ax.axhline(xmin = 0, xmax = t_dates[-1], ls = "--", linewidth = 2, color = "black")
    ax.axhline(xmin = 0, xmax = len(t_dates), ls = "--", linewidth = 2, color = "black")
    
    try:
        x_min = list(t_dates).index(str(sys.argv[5]))
        x_max = list(t_dates).index(str(sys.argv[6]))
    except:
        x_min = None
        x_max = None
    
    if (x_min is not None):
        ax.set_xlim((x_min, x_max))
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
    
    if x_min is not None:
        perday_orig = []
        for i in range(len(date_ticks)):
            perday_orig.append(list(t_dates).index(date_ticks[i]))
    else:
        perday_orig = perday
        
    ax.set_xticks(perday_orig)
    ax.set_xticklabels(date_ticks,
        rotation = 45, horizontalalignment = "right")

    #ax_twin.set_ylim((-0.02, 0.02))
    
    ax.legend(loc = (1.1, 0.75) ,fontsize = 20)
    pdf = PdfPages(sys.argv[w_save]+"/relative_fitness_groups.pdf")
    pdf.savefig(fig, bbox_inches = "tight")
    pdf.close()
 
    fig.savefig(sys.argv[w_save]+"/relative_fitness_groups.svg", bbox_inches = "tight")
    return status_list

num_groups = int(sys.argv[7])
w_save = 8
k = 9
ES_df_list = []
lineage_list = []
color_list = []
for i in range(num_groups):
    ES_df_list.append(pd.read_csv(ES_lin_dir+"/Susceptible_SpikeGroup_%s_all_PK.csv"%str(sys.argv[k+i])))    
    lineage_list.append(str(sys.argv[k+i]))
    color_list.append(str(sys.argv[k+num_groups+i]))

status_list = plot_fit(ES_df_list, lineage_list, color_list, w_save)
    
status = pd.DataFrame({"lineage":lineage_list, "relative_advantage":status_list})
status.to_csv(sys.argv[w_save]+"/plot_status.csv")
 

        