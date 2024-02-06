#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:28:29 2023

@author: raharinirina
"""

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
import pickle
import pdb
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
    ax.set_ylim((-0.1, 1.1*np.max(Z)))
    
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
PreFig(xsize = 20, ysize = 20)
fig = plt.figure(figsize = (10, 7))
ax = fig.add_subplot(1, 1, 1)

PK_df = pd.read_csv(sys.argv[1])
PK_df.drop(columns = "Unnamed: 0", inplace = True)
pk_cols = PK_df.columns
pk_t = PK_df["Day since activation"]

PK_all = PK_df.to_numpy()[:, pk_cols!="Day since activation"].astype(float) # all t_half, t_max simulations on the columns
PK_min = np.min(PK_all, axis = 1)
PK_max = np.max(PK_all, axis = 1)

Res_dir = sys.argv[3]
is_log=False
xval = "Days since antigen exposure"
yval = "Virus neutralization\n probability"
fig, ax = Display_Envelope(pk_t, np.array([PK_min]), np.array([PK_max]), 
                          is_log, labels = ["Epitopes PK (ranges)"], 
                          save_to = Res_dir+"/PK_Epitopes_ranges",
                          xval = xval, yval = yval,
                          linewidth = 4,
                          palette = "Greens",
                          alpha = 0.3,
                          fig = fig,
                          ax = ax,
                          antigen = "irrelevant")


fig.savefig(Res_dir+"/PK_Epitopes_ranges.svg", bbox_inches = "tight")

### Load P_Neut Data already pre-made ####
P_neut_dir = sys.argv[2]
num_groups = int(sys.argv[4])
k = 5
antigen = str(sys.argv[len(sys.argv)-1])

status = []
Lin_list = []
col_list = []
Lin_status = []

Min_list = []
Max_list =[]
xval = "Days since antigen exposure"
yval = "Virus neutralization\n probability"
s = 0
custom_col = sns.color_palette("Set2", 100) 

file = open("Spikegroups_membership.pck", "rb")
Pseudogroup_dic = pickle.load(file)
file.close()  
for i in range(num_groups):
    Lin_i_list = str(sys.argv[k+i])
    splited_var = np.array(Lin_i_list.split("/"))
    splited_var = splited_var[~(splited_var == "")]
    splited_var = splited_var[~(splited_var == " ")]
    for Lin_i in splited_var:
        try:
            try:
                Pneut_df = pd.read_csv(P_neut_dir+"/P_neut_%s.csv"%Lin_i)
                run=True
            except:
                try:
                    Pneut_df = pd.read_csv("results/Immunological_Landscape_ALL/P_neut_%s.csv"%Pseudogroup_dic[Lin_i])
                    run=True
                except:
                    try:
                        Pneut_df = pd.read_csv("results/Immunological_Landscape/P_neut_%s.csv"%Lin_i)
                        run = True
                    except:
                        print("Computation needed: P_neut file is not available for %s"%Lin_i)
                        run = False
            if run:
                try:
                    Pneut_df.drop(columns = "Unnamed: 0", inplace = True)
                except:
                    pass
                
                t = Pneut_df["Day since infection"] ### must be the same in all the Pneut files (which is the case in our pipeline)
    
                """Compute PNeut Envelope"""    
                EnvO_Min, EnvO_Max = Pneut_df["Proba Neut Min\n vs. %s antigen"%antigen].to_numpy(), Pneut_df["Proba Neut Max\n vs. %s antigen"%antigen].to_numpy()
                Min_list.append(EnvO_Min)
                Max_list.append(EnvO_Max)
                col_o = str(sys.argv[k+num_groups+i])
                if col_o in col_list:
                    if s<len(custom_col):
                        col_o = custom_col[s]
                    else:
                        rand_num = np.random.choice(1, 100)
                        col_o = sns.color_palette("rocked", rand_num)[0] 
                    s +=1
                    
                
                col_list.append(col_o)
                is_log=False
                
                ### save individual plots
                PreFig(xsize = 20, ysize = 20)
                figsize = (10,7)
                fig0 = plt.figure(figsize = figsize)
                ax0 = fig0.add_subplot(1, 1, 1)
                 
                fig0, ax0 = Display_Envelope(t, np.array([EnvO_Min]), np.array([EnvO_Max]), 
                                          is_log, 
                                          save_to = Res_dir+"/P_Neut_%s"%Lin_i,
                                          xval = xval, yval = yval,
                                          linewidth = 4,
                                          color = [col_o],
                                          alpha = 0.3,
                                          fig = fig0,
                                          ax = ax0,
                                          labels = [Lin_i],
                                          antigen = antigen) 
    
                status.append("Done")
                Lin_list.append(Lin_i)
            else:
                status.append("Not avail")
        except:
            status.append("Not avail")
            
        Lin_status.append(Lin_i)


PreFig(xsize = 20, ysize = 20)
figsize = (10,7)
fig = plt.figure(figsize = figsize)
ax = fig.add_subplot(1, 1, 1)
fig, ax = Display_Envelope(t, np.array(Min_list), np.array(Max_list), 
                          is_log, 
                          save_to = Res_dir+"/P_Neut_groups",
                          xval = xval, yval = yval,
                          linewidth = 4,
                          color = col_list,
                          alpha = 0.3,
                          fig = fig,
                          ax = ax,
                          labels = Lin_list,
                          antigen = antigen)  

pd.DataFrame({"Lineages":Lin_status, "status":status}).to_csv(Res_dir+"/plot_status.csv")
