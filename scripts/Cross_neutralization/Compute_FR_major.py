#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 22:07:05 2023

@author: raharinirina
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import numpy.ma as ma
import joblib as jb
from functools import partial
import re
import pickle
import sys
import pdb


"""Load SpikeGroups list"""
file1 = open(sys.argv[1], "rb") 
SpikeGroups_list = pickle.load(file1)["names"]
file1.close()

"""Load Mutation profile dictionary"""
file1 = open(sys.argv[2], "rb") 
mut_x_sites_dic = pickle.load(file1)["positions"]
file1.close()

variant_x_names_cross = SpikeGroups_list
# include wild-type
if "Wuhan-Hu-1" not in variant_x_names_cross:
    variant_x_names_cross = ["Wuhan-Hu-1"]+list(variant_x_names_cross)

mut_x_sites_dic["Wuhan-Hu-1"] = []

"""Load DMS Escape fraction data"""
Escape_Fraction = pd.read_csv(sys.argv[3])
Ab_classes = np.unique(Escape_Fraction["group"].values.astype(str))

"""Load lineage name to assess and it's mutation profile"""
try:
    Lin_name = sys.argv[4]
    mut_file = open(sys.argv[5], "r")
    mut_lin0 = mut_file.readlines()
    mut_file.close()

    mut_Lin = []
    for mut in mut_lin0:
        if mut[:3] not in ("DEL", "del"):
            if len(re.findall(r'\d+', mut))>0:
                mut_Lin.append(re.findall(r'\d+', mut)[0])       
                mut_Lin = list(np.unique(np.array(mut_Lin).astype(int)))

    """Update mutation profile dictionary"""
    mut_x_sites_dic_updated = mut_x_sites_dic.copy()
    mut_x_sites_dic_updated[Lin_name] = mut_Lin
except:
    mut_x_sites_dic_updated = mut_x_sites_dic.copy()

def sub_Bind(d, tiled_esc, Where_Mut, Where_Cond):
    Inter_Cond_Mut = Where_Mut & Where_Cond[np.newaxis, d, :]
    Bind_list_d = np.prod(1 - tiled_esc[:, np.newaxis, :]*Inter_Cond_Mut, axis = (1,2))
    Missing_cond_data_d = ~np.any(np.any(Where_Mut & Where_Cond[np.newaxis, d, :], axis = 2), axis = 1)
    return [Bind_list_d, Missing_cond_data_d]


def FR_xy(i, mut_sites, mut_bool_g1, mut_bool_g2, escape_ab_dic, ab, variant_name, EF_func = "MEAN", GM = False, quiet = True, joblib = None):
    vars_num = mut_bool_g2.shape[0]
    
    test_i = np.tile(mut_bool_g1[i, :], (vars_num, 1))
    diff_sites = (test_i ^ mut_bool_g2)

    escape_data_ab = escape_ab_dic["escape_data_ab"]
    
    conditions = escape_ab_dic["conditions"]
    ab_sub_list = escape_ab_dic["ab_sub_list"]   
    escape_sites = escape_ab_dic["escape_sites"]
    IC50_list = escape_ab_dic["IC50_list"]
    
    tiled_esc = np.tile(escape_data_ab, (vars_num, 1))
    Bind_list  = np.ones((vars_num, len(conditions)))
    Missing_cond_data = np.zeros((vars_num, len(conditions)), dtype = bool)
    Where_Cond = conditions[:, np.newaxis] == ab_sub_list[np.newaxis, :]
    tiled_mut = ma.array(np.tile(mut_sites, (vars_num, 1)), mask = ~diff_sites)
    Where_Mut = tiled_mut[:, :, np.newaxis] == escape_sites[np.newaxis, np.newaxis,  :]
    
    if joblib is not None:
        """Parallel codes --- macOS Monterey 12.5 crashes --- Not used by default """
        pfunc = partial(sub_Bind, tiled_esc = tiled_esc, Where_Mut = Where_Mut, Where_Cond = Where_Cond)
        try:
            jb_res = list(jb.Parallel(n_jobs = -1)(jb.delayed(pfunc)(d) for d in range(len(conditions))))
        except:
            jb_res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
        
        for d in range(len(conditions)):
            Bind_list[:, d]   = np.array(jb_res[d][0])
            Missing_cond_data[:, d] = np.array(jb_res[d][1])
    else:
        """ Brute force method """
        for d in range(len(conditions)):
            #print(d+1, len(conditions))
            Inter_Cond_Mut = Where_Mut & Where_Cond[np.newaxis, d, :]
            Bind_list[:, d] = np.prod(1 - tiled_esc[:, np.newaxis, :]*Inter_Cond_Mut, axis = (1,2))  
            Missing_cond_data[:, d] = ~np.any(np.any(Where_Mut & Where_Cond[np.newaxis, d, :], axis = 2), axis = 1)

    
    retained_binding = Bind_list
    FR_list = np.maximum((0.4/IC50_list[np.newaxis, :])*((1/retained_binding) - 1), 1)# DMS data was performed at fixed antibody concentration of 0.4 micro-gram/L (IC50 values are in micro-gram/L)
    #FR could be less than 1 happen if the IC50 is close to 0.4 and thus the binding retained is not reliable, we set those to 1
    FR_list = ma.array(FR_list, mask = Missing_cond_data) # do not take into account the missing data
    FR_ab = FR_list
    
    Missed = Missing_cond_data
    Greater_one = []
    
    return FR_ab, Missed, Greater_one 

def cross_reactivity(variant_name, escape_per_sites, Ab_classes, mut_sites_per_variant, EF_func = "MEAN", GM = False, quiet = True, joblib = None):
    FRxy = {}
    Missed = []
    Greater_one = []
    ### extract variants groups
    variants_g1, variants_g2 = variant_name
    ### extract all sites
    mut_sites = []
    for var in list(mut_sites_per_variant.keys()):
        if (var in variants_g1) or (var in variants_g2):
            mut_sites += list(np.array(mut_sites_per_variant[var]).astype(str))
    
    mut_sites = np.unique(mut_sites).astype(str)
    ### Construct a boolean array for location of mutation for all variants    
    mut_bool_g1 = np.zeros((len(variants_g1), len(mut_sites)), dtype = bool)
    mut_bool_g2= np.zeros((len(variants_g2), len(mut_sites)), dtype = bool)
    #print(variant_name)
    for i in range(max(len(variants_g1), len(variants_g2))):
        for j in range(len(mut_sites)): 
            if i<len(variants_g1):
                mut_bool_g1[i, j] = mut_sites[j] in list(np.array(mut_sites_per_variant[variants_g1[i]]).astype(str))
            if i<len(variants_g2):
                mut_bool_g2[i, j] = mut_sites[j] in list(np.array(mut_sites_per_variant[variants_g2[i]]).astype(str))

    for ab in Ab_classes:
        FRxy_ab = np.ones((len(variants_g1), len(variants_g2)))
        Missing_Cond = np.zeros((len(variants_g1), len(variants_g2)))
        
        escape_ab_dic = {}
        escape_data = escape_per_sites["mut_escape"].values
        where_ab_group = (escape_per_sites["group"].values).astype(str) == str(ab)
        escape_ab_dic["ab_sub_list"]   = (escape_per_sites["condition"].values).astype(str)[where_ab_group]

        sub_table = escape_per_sites.loc[(escape_per_sites['condition'][where_ab_group].drop_duplicates().index), ['condition', 'IC50']]
        escape_ab_dic["conditions"] = sub_table["condition"].values
        escape_ab_dic["IC50_list"] = sub_table["IC50"].values
        
        escape_ab_dic["escape_sites"] = ((escape_per_sites["site"].values).astype(str))[where_ab_group]
        
        escape_ab_dic["escape_data_ab"] = escape_data[where_ab_group]
        
        for i in range(len(variants_g1)):
            FR, missed, gOne = FR_xy(i, mut_sites, mut_bool_g1, mut_bool_g2, escape_ab_dic, ab, variant_name, mut_sites_per_variant, GM = GM, quiet = quiet, joblib = joblib)
            if EF_func == "MEAN":
                FRxy_ab[i, :] = np.mean(FR, axis = 1)
                
            elif EF_func == "PROD":
                FRxy_ab[i, :] = np.prod(FR, axis = 1)
            
            Missing_Cond[i, :] = np.all(missed, axis = 1)
            
            Missed = []
            Greater_one += gOne
            
        #FRxy[ab] = ma.array(FRxy_ab, mask = Missing_Cond.astype(bool), fill_value = np.nan)
        FRxy_ab[Missing_Cond.astype(bool)] = 1
        FRxy[ab] = FRxy_ab
        
    Missed = np.unique(np.array(Missed))
    Greater_one = np.unique(np.array(Greater_one))           
    return FRxy, Missed, Greater_one

 
### Compute Cross reactivity between major variant groups for sanity checks, 
#only computed when the timeline is wide enough to contain the major variant groups
file = open("Spikegroups_membership.pck", "rb")
Pseudogroup_dic = pickle.load(file)
file.close()

k = 7
lineages_sim = []
while k<(int(sys.argv[6])+7):
    lineages_sim.append(sys.argv[k])
    k +=1
    
mut_sim = []
while k<(int(sys.argv[6])+7+len(lineages_sim)):
    mut_sim.append(sys.argv[k])
    k +=1
        
Top_Pseudo = []
Pseudo_keys = list(Pseudogroup_dic.keys())

for spklin in lineages_sim:
    if spklin in Pseudo_keys:
        Top_Pseudo.append(spklin)
        mut_maj = mut_x_sites_dic[Pseudogroup_dic[spklin]]
        """Update mutation profile dictionary"""        
        mut_x_sites_dic_updated[spklin] = mut_maj
    else:           
        try:
            mut_file = open(mut_sim[lineages_sim.index(spklin)], "r")
            mut_lin0 = mut_file.readlines()
            mut_file.close()
            mut_maj= []
            for mut in mut_lin0:
                if mut[:3] not in ("DEL", "del"):
                    if len(re.findall(r'\d+', mut))>0:
                        mut_maj.append(re.findall(r'\d+', mut)[0])       
                        mut_maj = list(np.unique(np.array(mut_maj).astype(str)))
            
            """Update mutation profile dictionary"""        
            mut_x_sites_dic_updated[spklin] = mut_maj
            Top_Pseudo.append(spklin)
        except:
            pass

if "Wuhan-Hu-1" not in Top_Pseudo:
    Top_Pseudo = ["Wuhan-Hu-1"] + list(Top_Pseudo)
    
a = 1
if len(Top_Pseudo)!=0:
    Cross_react_dic = {}
    mut_x_sites_dic_used = mut_x_sites_dic_updated.copy()
    try:
        Top_Pseudo.append(Lin_name)
        mut_x_sites_dic_used[Lin_name] = mut_Lin
    except:
        pass
      
    for ab in Ab_classes:
        print("Assess major lineages/pseudogroups with the NTD-RBD mutation positions ")
        print("Cross reactivity countdown", a, "out of %d epitope clases"%len(Ab_classes))
        if ab!= "NTD":
            Cross_Lin, Missed, Greater_one = cross_reactivity((Top_Pseudo, Top_Pseudo), 
                       Escape_Fraction, 
                       [ab],
                       mut_x_sites_dic_used)
            
            FRxy_ab = Cross_Lin[ab]
            Cross_react_dic[ab] = FRxy_ab
            
        a +=1
    
    Cross_react_dic["variant_list"] = list(Top_Pseudo)
    
    """Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""   
    n = len(Cross_react_dic["variant_list"])
    FR_NTB = np.ones((n, n))
    for i in range(n):
        var_1 = Cross_react_dic["variant_list"][i]
        for j in range(n):
            if i > j:
                var_2 = Cross_react_dic["variant_list"][j]
    
                sites_1 = set(np.array(mut_x_sites_dic_used[var_1]).astype(int))
                sites_2 = set(np.array(mut_x_sites_dic_used[var_2]).astype(int))
    
                sites = list(sites_1.symmetric_difference(sites_2))
                FR_sites = 1
                for s in sites:
                    s = int(s)
                    if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                        FR_sites *= 10
                FR_NTB[i, j] = FR_sites
                FR_NTB[j, i] = FR_sites
        
    Cross_react_dic["NTD"] = FR_NTB
    file0 = open(sys.argv[k], "wb") 
    pickle.dump(Cross_react_dic, file0)
    file0.close()    
