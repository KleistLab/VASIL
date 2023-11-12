#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 10:03:18 2023

@author: raharinirina
"""
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
#### Visualisation ###  
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)


def Display_Envelope(t, Y, Z, is_log, labels, figsize = (7, 7), xysize = (15,15), labsize = 20, save_to = "test", xval = "x", yval = "f(x)", 
                    linewidth = 3, palette = None, linestyle = None, color = None, alpha = 0.5, ax = None, fig = None, antigen = "(Not specified)"):
    if fig == None:
        PreFig(xsize = xysize[0], ysize = xysize[1])
        fig = plt.figure(figsize = figsize)
        if ax == None:
            ax = fig.add_subplot(1, 1, 1)
    
    if linestyle is None:
        linestyle = ["-"]*Y.shape[0]
    
    if palette is None:
        for i in range(Y.shape[0]):
            if color is None:
                ax.fill_between(t, Y[i, :], Z[i, :], label = labels[i], alpha = alpha)
            else:
                ax.fill_between(t, Y[i, :], Z[i, :], label = labels[i], alpha = alpha, color = color[i])
    else:
        col = sns.color_palette(palette, Y.shape[0])
        for i in range(Y.shape[0]):
            ax.fill_between(t, Y[i, :], Z[i, :], label = labels[i], alpha = alpha, color = col[i])
        
    if is_log:
        plt.ylabel("$\ln$ %s"%yval, fontsize = labsize)
    else:
        plt.ylabel("%s"%yval, fontsize = labsize)  
    
    plt.xlabel(xval, fontsize = labsize)
    
    if labels != [""]:
        ax.legend(loc = (1.2, 0.) ,fontsize = labsize, ncols = np.ceil(len(labels)/4))
    
    if save_to is not None:
        ### save figure in pdf ###
        if is_log:
             ax.set_ylabel("$\ln$ %s"%yval, fontsize = labsize)
        else:
             ax.set_ylabel("%s"%yval, fontsize = labsize)
             
        ax.set_xlabel(xval, fontsize = labsize)
        if antigen != "irrelevant":
            ax.set_title("%s antigen"%antigen, fontsize = labsize)
            
        pdf = PdfPages(save_to+".pdf")
        pdf.savefig(fig, bbox_inches = "tight")
        pdf.close()
        
        ### save figure as svg
        fig.savefig(save_to+".svg", bbox_inches = "tight")
    
    return fig, ax

### Load PK data already pre-made ###
PK_df = pd.read_csv(sys.argv[1])
PK_df.drop(columns = "Unnamed: 0", inplace = True)
pk_cols = PK_df.columns
pk_t = PK_df["Day since activation"]

PK_all = PK_df.to_numpy()[:, pk_cols!="Day since activation"].astype(float) # all t_half, t_max simulations on the columns
PK_min = np.min(PK_all, axis = 1)
PK_max = np.max(PK_all, axis = 1)

PreFig(xsize = 20, ysize = 20)
fig = plt.figure(figsize = (10, 7))
ax = fig.add_subplot(1, 1, 1)
is_log=False
xval = "Days since antigen exposure"
yval = "Virus neutralization\n probability"
fig, ax = Display_Envelope(pk_t, np.array([PK_min]), np.array([PK_max]), 
                          is_log, labels=["Epitopes PK (ranges)"], 
                          save_to = sys.argv[5]+"/PK_Epitopes_ranges",
                          xval = xval, yval = yval,
                          linewidth = 4,
                          palette = "Greens",
                          alpha = 0.3,
                          fig = fig,
                          ax = ax,
                          antigen = "irrelevant")

fig.savefig(sys.argv[5]+"/PK_Epitopes_ranges.svg", bbox_inches = "tight")
status = {"PK":["done"]}

### Load P_Neut Data already pre-made ####
Pneut_df = pd.read_csv(sys.argv[2])
Pneut_df.drop(columns = "Unnamed: 0", inplace = True)
t = Pneut_df["Day since infection"]

PreFig(xsize = 20, ysize = 20)
figsize = (10,7)
fig = plt.figure(figsize = figsize)
ax = fig.add_subplot(1, 1, 1)

"""Compute PNeut Envelope"""
antigen = str(sys.argv[6])    
EnvO_Min, EnvO_Max = Pneut_df["Proba Neut Min\n vs. %s antigen"%antigen].to_numpy(), Pneut_df["Proba Neut Max\n vs. %s antigen"%antigen].to_numpy()
col_o = str(sys.argv[3])
is_log=False
xval = "Days since antigen exposure"
yval = "Virus neutralization\n probability"
Lin = str(sys.argv[4])
fig, ax = Display_Envelope(t, np.array([EnvO_Min]), np.array([EnvO_Max]), 
                          is_log, labels = [Lin],
                          save_to = sys.argv[5]+"/P_Neut_"+Lin,
                          xval = xval, yval = yval,
                          linewidth = 4,
                          color = [col_o],
                          alpha = 0.3,
                          fig = fig,
                          ax = ax,
                          antigen = antigen) 

status["P_neut"] = ["done"]
fig.savefig(sys.argv[5]+"/P_Neut_"+str(sys.argv[4])+".svg", bbox_inches = "tight")
(pd.DataFrame(status)).to_csv(sys.argv[5]+"/plot_status.csv")
