#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pickle
import sys
import seaborn as sns
import numpy.ma as ma
import pdb


#### Load Data #####
cross_file = open(sys.argv[1], "rb")
Cross_Epitope_Dic_orig = pickle.load(cross_file)
cross_file.close()
file = open("Spikegroups_membership.pck", "rb")
Pseudogroup_dic = pickle.load(file)
file.close()

try:
    lineages_sim = ["BA.2", "BA.4", "BA.5", "BQ.1.1", "BE.1.1", "CH.1.1", "XBB.1.5"]
    Top_Pseudo = []
    Top_lab = []
    Pseudo_done = []
    for spklin in lineages_sim:
        if Pseudogroup_dic[spklin] not in Pseudo_done:
            Top_Pseudo.append(Pseudogroup_dic[spklin])
            Pseudo_done.append(Pseudogroup_dic[spklin])
            if Pseudogroup_dic[spklin] != spklin:
                Top_lab.append(Pseudogroup_dic[spklin]+"/"+spklin)
            else:
                Top_lab.append(spklin)
        else:
            ix = Pseudo_done.index(Pseudogroup_dic[spklin])
            Top_lab[ix] = Top_lab[ix]+"/"+spklin
           
    
    Top_Pseudo = ["Wuhan-Hu-1"] + list(Top_Pseudo)
    Top_lab = ["Wuhan-Hu-1"]+Top_lab
    Pseudo_lab_cross = Top_lab

except:
    Top_Pseudo = np.array(Cross_Epitope_Dic_orig["variant_list"][:10])
    if "Wuhan-Hu-1" not in Top_Pseudo:
        Top_Pseudo = ["Wuhan-Hu-1"] + list(Top_Pseudo[:-1])
    else:
        Top_Pseudo = ["Wuhan-Hu-1"] + list(Top_Pseudo[Top_Pseudo!="Wuhan-Hu-1"])
        
    Pseudo_lab_cross = []
    for i in range(len(Top_Pseudo)):
        if Top_Pseudo[:7] == "Spike. ":
            Pseudo_lab_cross.append(Top_Pseudo[i][7:])
        else:
            Pseudo_lab_cross.append(Top_Pseudo[i])


All_Pseudo = list(Cross_Epitope_Dic_orig["variant_list"])
Cross_Epitope_Dic_orig.pop("variant_list")
choosen_Ab = list(Cross_Epitope_Dic_orig.keys())

Cross_Epitope_Dic = {}
for ab in choosen_Ab:
    Cross = np.ones((len(Top_Pseudo),len(Top_Pseudo)))
    for i in range(len(Top_Pseudo)):
        w_i = All_Pseudo.index(Top_Pseudo[i])
        for j in range(len(Top_Pseudo)):
            w_j = All_Pseudo.index(Top_Pseudo[j])
            Cross[i,j] = Cross_Epitope_Dic_orig[ab][w_i, w_j]
    Cross_Epitope_Dic[ab] = Cross
        
### Colors picked from figure sketch.pptx 3D structure of epitope classes ### Replace RGBA ####
#cmap_base = ["Reds", "Blues", "Wistia", "Wistia"]
from matplotlib.colors import ListedColormap
Reds_0 = [(255/255., 248/255., 247/255.) , (253/255., 242/255., 242/255.), (250/255., 216/255., 213/255.), 
          (245/255., 175/255., 169/255.), (237/255., 105/255., 95/255.) , (231/255., 54/255., 40/255.)] # for normal scale

Reds_1 = [(255/255., 248/255., 247/255.), (250/255., 216/255., 213/255.), 
          (245/255., 175/255., 169/255.), (231/255., 54/255., 40/255.)] # for log10 scale

Reds = Reds_1
cmReds = ListedColormap(Reds, N = len(Reds)) 

Blues_1 = [(250/255., 250./255, 1.), (232/255., 231/255., 1.), 
           (134/255., 130./255., 253/255.), (37/255., 26/255, 252/255.)] # for log10 scale
Blues = Blues_1
cmBlues = ListedColormap(Blues, N= len(Blues))

Yellows_1 = [(255/255., 255./255, 231/255.), (254/255., 254/255., 190/255.),
            (252/255., 252/255, 6./255.),
             (252/255., 237/255., 100./255), (252/255, 219/255., 96/255.)] # for log10 scale

Yellows = Yellows_1
cmYellows = ListedColormap(Yellows, N= len(Yellows))

cmap_dic = {"A":cmReds, "B":cmBlues, "D1":cmYellows, "D2":cmYellows}


cbar_labsize = 50
title_labsize = 50
ticksize = (50, 50)
facecolor = "white"

#### Visualisation ###  
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)


xysize = (15, 15) # size of xy ticks

try:
    Lin_name = sys.argv[2]
    if Lin_name not in Top_Pseudo:
        cross_lin = open(sys.argv[3], "rb")
        Cross_Epitope_lin = pickle.load(cross_lin)
        var_lin = Cross_Epitope_lin["variant_list"]
        cross_lin.close()
        for ab in choosen_Ab:
            FR_lin = np.ones(len(Top_Pseudo)+1)
            for i in range(len(Top_Pseudo)):
                FR_lin[1+i] = Cross_Epitope_lin[ab][list(var_lin).index(Lin_name), list(var_lin).index(Top_Pseudo[i])]
                
            Cross_Epitope_Dic[ab] = np.row_stack((Cross_Epitope_Dic[ab], FR_lin[1:]))
            FR_linT = np.ones(len(Top_Pseudo)+1)
            FR_linT[:-1] = FR_lin[1:]
            Cross_Epitope_Dic[ab] = np.column_stack((Cross_Epitope_Dic[ab], FR_linT.T))
        
        Top_Pseudo = list(Top_Pseudo) + [Lin_name] 
        Pseudo_lab_cross = list(Pseudo_lab_cross) + [Lin_name]
except:
    try:
        Lin_name = sys.argv[2]
        cross_lin = open(sys.argv[3], "rb")
        Cross_Epitope_lin = pickle.load(cross_lin)
        var_lin = Cross_Epitope_lin["variant_list"]
        Top_Pseudo = var_lin[:9]
        
        Cross_Epitope_Dic = {}
        if "Wuhan-Hu-1" not in Top_Pseudo:
            Top_Pseudo = ["Wuhan-Hu-1"] + list(Top_Pseudo[:-1])
        else:
            Top_Pseudo = ["Wuhan-Hu-1"] + list(Top_Pseudo[Top_Pseudo!="Wuhan-Hu-1"])
        
        Pseudo_lab_cross = []
        for i in range(len(Top_Pseudo)):
            if Top_Pseudo[:7] == "Spike. ":
                Pseudo_lab_cross.append(Top_Pseudo[i][7:])
            else:
                Pseudo_lab_cross.append(Top_Pseudo[i])
                
        if Lin_name not in Top_Pseudo:            
            Top_Pseudo = list(Top_Pseudo) + [Lin_name] 
            Pseudo_lab_cross = list(Pseudo_lab_cross) + [Lin_name]
            
        for ab in choosen_Ab:
            Cross = np.ones((len(Top_Pseudo), len(Top_Pseudo)))
            for i in range(len(Top_Pseudo)):
                for j in range(len(Top_Pseudo)):
                    if (Top_Pseudo[i] != Lin_name)&(Top_Pseudo[j]!= Lin_name):
                        w_i = list(All_Pseudo).index(Top_Pseudo[i])
                        w_j = list(All_Pseudo).index(Top_Pseudo[j])
                        Cross[i, j] = Cross_Epitope_Dic_orig[ab][w_i, w_j]
                    else:
                        w_i = list(var_lin).index(Top_Pseudo[i])
                        w_j = list(var_lin).index(Top_Pseudo[j])
                        Cross[i, j] = Cross_Epitope_lin[ab][w_i, w_j]
            Cross_Epitope_Dic[ab] = Cross
                
    except:
        pass


triup = np.triu_indices(len(Top_Pseudo), k=0)
mask_triup = np.ones((len(Top_Pseudo), len(Top_Pseudo))).astype(bool)
mask_triup[triup] = False
mask_triup = mask_triup.astype(bool)

for k in range(len(choosen_Ab)):
    ab = choosen_Ab[k]
    Cross_ab = Cross_Epitope_Dic[ab]
    ### Set Colorbar limit ###
    Cross_sub = np.log10(Cross_ab)
    
    if ab == "A":
        FR_vals_sub = np.arange(0, np.max(Cross_sub)+0.5, 0.5)
        center = 1
    elif ab == "B":
        FR_vals_sub = np.arange(0, np.max(Cross_sub)+1, 1)
        center = 1.*np.mean(FR_vals_sub)
    elif ab == "D1":
        FR_vals_sub = np.arange(0, np.max(Cross_sub)+0.5, 0.5)
        center = 1.*np.mean(FR_vals_sub)
    else:
        FR_vals_sub = np.arange(0, np.max(Cross_sub)+1, 1)
        center = 1.*np.mean(FR_vals_sub)
    
    ### show >= ... only if there is enough distance between max and the last FR_vals (otherwise it is crowded)
    if FR_vals_sub[-1] < np.max(Cross_sub) - 0.1:
        FR_vals_sub = np.append(FR_vals_sub[:-1], np.max(Cross_sub))
  
    #Cross_sub[Cross_sub>=FR_vals_sub[-1]] = FR_vals_sub[-1]
    
    ### Put mask only after setting colorbar limits ####
    Cross_sub = ma.array(Cross_sub, mask = mask_triup)

    Cross_Epitope_Dic["Resistance to %s"%ab] = Cross_sub
    
    PreFig(xsize = xysize[0], ysize = xysize[1])
    fig_fr = plt.figure(figsize = (16, 13))
    ax_fr = fig_fr.add_subplot(1, 1, 1)
    
    dLab = "Cross-Resistance to %s"%ab
    
    if ab in list(cmap_dic.keys()):
        cMap = sns.heatmap(data = Cross_sub,
                    mask = mask_triup, 
                    cmap = cmap_dic[ab], 
                    xticklabels = Pseudo_lab_cross, 
                    yticklabels = Pseudo_lab_cross, 
                    cbar = True, 
                    annot = True, 
                    #annot = False,
                    annot_kws = {"size": 30},
                    fmt = ".2f", ## annotations decimals
                    cbar_kws = {'label': 'FR (log 10)',"shrink": 0.75, "ticks":FR_vals_sub}, 
                    center = center,
                    linewidths = 3,
                    linecolor = "white")
    else:
        cMap = sns.heatmap(data = Cross_sub,
                    mask = mask_triup, 
                    xticklabels = Pseudo_lab_cross, 
                    yticklabels = Pseudo_lab_cross, 
                    cbar = True, 
                    annot = True, 
                    #annot = False,
                    annot_kws = {"size": 30},
                    fmt = ".2f", ## annotations decimals
                    cbar_kws = {'label': 'FR (log 10)',"shrink": 0.75, "ticks":FR_vals_sub}, 
                    center = center,
                    linewidths = 3,
                    linecolor = "white")
        
    cbar = cMap.figure.axes[-1] # get colorbar instance
    cbar.yaxis.label.set_size(cbar_labsize)
                
    cbar.tick_params(labelsize = cbar_labsize)
    cMap.set_facecolor(facecolor) 
    plt.title(dLab, fontsize = title_labsize)
    ax_fr.set_aspect("equal")
    plt.xticks(fontsize = ticksize[0], rotation = 45, horizontalalignment = "right")
    plt.yticks(fontsize = ticksize[1], rotation = 0)
    
    plt.subplots_adjust(hspace=0.65, wspace=0.5)
    
    cbar = ax_fr.collections[0].colorbar
    cbar.set_ticklabels(["%.1f"%FR_vals_sub[:-1][i] for i in range(len(FR_vals_sub[:-1]))]+["$\geq$"+" %.1f"%FR_vals_sub[-1]])
    
    pdf = PdfPages(sys.argv[4]+"/Cross_React_AB_%s.pdf"%ab)
    pdf.savefig(fig_fr, bbox_inches = "tight")
    fig_fr.savefig(sys.argv[4]+"/Cross_React_AB_%s.svg"%ab, bbox_inches = "tight")
    pdf.close()
 
### Always plot Cross reactivity between major variant groups for sanity checks, 
#only computed when the timeline is wide enough to contain the major variant groups   
try:    
    file0 = open("results/Cross_to_major_variants.pck", "rb") 
    Cross_show=pickle.load(file0)
    file0.close()
    Top_Pseudo_var = Cross_show["variant_list"]
    Top_Pseudo = []
    Top_lab = []
    Pseudo_done = []
    lineages_sim = lineages_sim + [Lin_name]
    for spklin in lineages_sim:
        if (spklin != Lin_name):
            if (Pseudogroup_dic[spklin] in Top_Pseudo_var):
                if Pseudogroup_dic[spklin] not in Pseudo_done:
                    Top_Pseudo.append(Pseudogroup_dic[spklin])
                    Pseudo_done.append(Pseudogroup_dic[spklin])
                    if Pseudogroup_dic[spklin] != spklin:
                        Top_lab.append(Pseudogroup_dic[spklin]+"/"+spklin)
                    else:
                        Top_lab.append(spklin)
                else:
                    ix = Pseudo_done.index(Pseudogroup_dic[spklin])
                    Top_lab[ix] = Top_lab[ix]+"/"+spklin
        else:
            Top_Pseudo.append(spklin)
            Top_lab.append(spklin)   
            
    if "Wuhan-Hu-1" not in Top_Pseudo:
        Top_Pseudo = ["Wuhan-Hu-1"] + list(Top_Pseudo)
        Top_lab = ["Wuhan-Hu-1"]+Top_lab
    Pseudo_lab_cross = Top_lab
    Cross_Dic = {}
    
    for ab in choosen_Ab:
        Cross = np.ones((len(Top_Pseudo),len(Top_Pseudo)))
        for i in range(len(Top_Pseudo)):
            w_i = Top_Pseudo_var.index(Top_Pseudo[i])
            for j in range(len(Top_Pseudo)):
                w_j = Top_Pseudo_var.index(Top_Pseudo[j])
                Cross[i,j] = Cross_show[ab][w_i, w_j]
        Cross_Dic[ab] = Cross
        
    triup = np.triu_indices(len(Top_Pseudo), k=0)
    mask_triup = np.ones((len(Top_Pseudo), len(Top_Pseudo))).astype(bool)
    mask_triup[triup] = False
    mask_triup = mask_triup.astype(bool)

    for k in range(len(choosen_Ab)):
        ab = choosen_Ab[k]
        Cross_ab = Cross_Dic[ab]
        ### Set Colorbar limit ###
        Cross_sub = np.log10(Cross_ab)
        
        if ab == "A":
            FR_vals_sub = np.arange(0, np.max(Cross_sub)+0.5, 0.5)
            center = 1
        elif ab == "B":
            FR_vals_sub = np.arange(0, np.max(Cross_sub)+1, 1)
            center = 1.*np.mean(FR_vals_sub)
        elif ab == "D1":
            FR_vals_sub = np.arange(0, np.max(Cross_sub)+0.5, 0.5)
            center = 1.*np.mean(FR_vals_sub)
        else:
            FR_vals_sub = np.arange(0, np.max(Cross_sub)+1, 1)
            center = 1.*np.mean(FR_vals_sub)
        
        ### show >= ... only if there is enough distance between max and the last FR_vals (otherwise it is crowded)
        if FR_vals_sub[-1] < np.max(Cross_sub) - 0.1:
            FR_vals_sub = np.append(FR_vals_sub[:-1], np.max(Cross_sub))
      
        #Cross_sub[Cross_sub>=FR_vals_sub[-1]] = FR_vals_sub[-1]
        
        ### Put mask only after setting colorbar limits ####
        Cross_sub = ma.array(Cross_sub, mask = mask_triup)

        Cross_Epitope_Dic["Resistance to %s"%ab] = Cross_sub
        
        PreFig(xsize = xysize[0], ysize = xysize[1])
        fig_fr = plt.figure(figsize = (16, 13))
        ax_fr = fig_fr.add_subplot(1, 1, 1)
        
        dLab = "Cross-Resistance to %s"%ab
        
        if ab in list(cmap_dic.keys()):
            cMap = sns.heatmap(data = Cross_sub,
                        mask = mask_triup, 
                        cmap = cmap_dic[ab], 
                        xticklabels = Pseudo_lab_cross, 
                        yticklabels = Pseudo_lab_cross, 
                        cbar = True, 
                        annot = True, 
                        #annot = False,
                        annot_kws = {"size": 30},
                        fmt = ".2f", ## annotations decimals
                        cbar_kws = {'label': 'FR (log 10)',"shrink": 0.75, "ticks":FR_vals_sub}, 
                        center = center,
                        linewidths = 3,
                        linecolor = "white")
        else:
            cMap = sns.heatmap(data = Cross_sub,
                        mask = mask_triup, 
                        xticklabels = Pseudo_lab_cross, 
                        yticklabels = Pseudo_lab_cross, 
                        cbar = True, 
                        annot = True, 
                        #annot = False,
                        annot_kws = {"size": 30},
                        fmt = ".2f", ## annotations decimals
                        cbar_kws = {'label': 'FR (log 10)',"shrink": 0.75, "ticks":FR_vals_sub}, 
                        center = center,
                        linewidths = 3,
                        linecolor = "white")
            
        cbar = cMap.figure.axes[-1] # get colorbar instance
        cbar.yaxis.label.set_size(cbar_labsize)
        cbar.tick_params(labelsize = cbar_labsize)
        cMap.set_facecolor(facecolor) 
        plt.title(dLab, fontsize = title_labsize)
        ax_fr.set_aspect("equal")
        plt.xticks(fontsize = ticksize[0], rotation = 45, horizontalalignment = "right")
        plt.yticks(fontsize = ticksize[1], rotation = 0)
        
        plt.subplots_adjust(hspace=0.65, wspace=0.5)
        
        cbar = ax_fr.collections[0].colorbar
        cbar.set_ticklabels(["%.1f"%FR_vals_sub[:-1][i] for i in range(len(FR_vals_sub[:-1]))]+["$\geq$"+" %.1f"%FR_vals_sub[-1]])
        
        pdf = PdfPages(sys.argv[4]+"/major_Cross_React_AB_%s.pdf"%ab)
        pdf.savefig(fig_fr, bbox_inches = "tight")
        fig_fr.savefig(sys.argv[4]+"/major_Cross_React_AB_%s.svg"%ab, bbox_inches = "tight")
        pdf.close()
except:
    pass

