#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pickle


from .util import cross_reactivity

### Load SpikeGroups list
file1 = open("Data/SpikeGroups.pck", "rb") 
SpikeGroups_list = pickle.load(file1)["names"]
file1.close()

### Load Mutation profile dictionary
file1 = open("Data/Mutation_Profiles.pck", "rb") 
mut_x_sites_dic = pickle.load(file1)["positions"]
file1.close()

### Load DMS data 
Escape_Fraction = pd.read_csv("Data/dms_per_ab_per_site.csv")
Ab_classes = np.unique(Escape_Fraction["group"])

start = 0
stop = len(SpikeGroups_list) - 1
variant_x_names_cross = SpikeGroups_list[start:stop+1]
# include wild-type
if "Wuhan-Hu-1" not in variant_x_names_cross:
    variant_x_names_cross = ["Wuhan-Hu-1"]+list(variant_x_names_cross)

mut_x_sites_dic["Wuhan-Hu-1"] = []

Cross_react_dic = {}
AB = ""
for ab in Ab_classes:
    Cross_react_dic_ab, Missed, Greater_one = cross_reactivity(variant_x_names_cross, Escape_Fraction, [ab], 
                                                          mut_x_sites_dic)
    
    Cross_react_dic_ab["variant_list"] = variant_x_names_cross
    Cross_react_dic[ab] = Cross_react_dic_ab[ab]

"""Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""   
n = len(variant_x_names_cross)
FR_NTB = np.ones((n, n))
for i in range(n):
    var_1 = variant_x_names_cross[i]
    for j in range(n):
        if i > j:
            var_2 = variant_x_names_cross[j]

            sites_1 = set(np.array(mut_x_sites_dic[var_1]).astype(int))
            sites_2 = set(np.array(mut_x_sites_dic[var_2]).astype(int))

            sites = list(sites_1.symmetric_difference(sites_2))
            FR_sites = 1
            for s in sites:
                s = int(s)
                if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                    FR_sites *= 10
            FR_NTB[i, j] = FR_sites
            FR_NTB[j, i] = FR_sites
        
Cross_react_dic["NTD"] = FR_NTB

Cross_react_dic["variant_list"] = variant_x_names_cross
file0 = open("Data/Cross_react_dic_spikegroups_%d_%d.pck"%(start, stop), "wb") 
pickle.dump(Cross_react_dic, file0)
file0.close()


### Cross_reactivity containing delta validation generated previously
file1 = open("Cross_react_dic_show.pck", "rb") # premade simulations
Cross_react_dic_show = pickle.load(file1)
variant_x_names_show = Cross_react_dic_show["variant_list"]
file1.close()

"""Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""   
n = len(variant_x_names_show)
FR_NTB = np.ones((n, n))
for i in range(n):
    var_1 = variant_x_names_show[i]
    for j in range(n):
        if i > j:
            var_2 = variant_x_names_show[j]

            sites_1 = set(np.array(mut_x_sites_dic[var_1]).astype(int))
            sites_2 = set(np.array(mut_x_sites_dic[var_2]).astype(int))

            sites = list(sites_1.symmetric_difference(sites_2))
            FR_sites = 1
            for s in sites:
                s = int(s)
                if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                    FR_sites *= 10
            FR_NTB[i, j] = FR_sites
            FR_NTB[j, i] = FR_sites

Cross_react_dic_show["NTD"] = FR_NTB
Cross_with_delta_validation = Cross_react_dic_show
Cross_react_dic_show["variant_list"] = variant_x_names_show
file0 = open("Data/Cross_with_delta_validation.pck", "wb") 
pickle.dump(Cross_with_delta_validation, file0)
file0.close()
