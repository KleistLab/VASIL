#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pickle
import sys

from .util import cross_reactivity

### Load SpikeGroups list
file1 = open(sys.argv[1], "rb") 
SpikeGroups_list = pickle.load(file1)["names"]
file1.close()

### Load Mutation profile dictionary
file1 = open(sys.argv[2], "rb") 
mut_x_sites_dic = pickle.load(file1)["positions"]
file1.close()

### Load DMS data 
Escape_Fraction = pd.read_csv(sys.argv[3])
Ab_classes = np.unique(Escape_Fraction["group"].values.astype(str))

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
file0 = open(sys.argv[4], "wb") 
pickle.dump(Cross_react_dic, file0)
file0.close()


### Compute Cross reactivity to Delta for valitation
# use Delta: B.1.617.2 for validation (Used in Clinical data paper)
variant_x_names_show = ["Wuhan-Hu-1", "Delta: B.1.617.2"]
mut_dic_show = {"Wuhan-Hu-1":[], "Delta: B.1.617.2": [614, 950, 142, 452, 681, 19, 478]}

Cross_with_delta_validation, Missed, Greater_one = cross_reactivity(variant_x_names_show, 
															 Escape_Fraction, Ab_classes[:-1], 
                                                        	 mut_dic_show)

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

Cross_with_delta_validation["NTD"] = FR_NTB
Cross_with_delta_validation["variant_list"] = variant_x_names_show
file0 = open(sys.argv[5], "wb") 
pickle.dump(Cross_with_delta_validation, file0)
file0.close()
