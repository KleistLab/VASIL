#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pickle


### Load SpikeGroups and Mutation profile 
from .util import cross_reactivity

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
    func_type = "MEAN"
    Cross_react_dic_ab, Missed, Greater_one = cross_reactivity(variant_x_names_cross, Escape_Fraction, [ab], 
                                                          mut_x_sites_dic, EF_func = func_type, quiet = False)
    
    Cross_react_dic_ab["variant_list"] = variant_x_names_cross
    Cross_react_dic_ab["func_type"]= func_type
    Cross_react_dic[ab] = Cross_react_dic_ab[ab]

Cross_react_dic["variant_list"] = variant_x_names_cross
Cross_react_dic["func_type"]= func_type
file0 = open("Data/Cross_react_dic_pseudogroups_%d_%d.pck"%(start, stop), "wb") # <--- check FILENAME
pickle.dump(Cross_react_dic, file0)
file0.close()


