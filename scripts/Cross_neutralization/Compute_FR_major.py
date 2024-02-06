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
Mut_infos = pickle.load(file1)
mut_x_sites_dic = Mut_infos["positions"]
AA_change_dic = Mut_infos["AA_changes"]
file1.close()

variant_x_names_cross = SpikeGroups_list
# include wild-type
if "Wuhan-Hu-1" not in variant_x_names_cross:
    variant_x_names_cross = ["Wuhan-Hu-1"]+list(variant_x_names_cross)

mut_x_sites_dic["Wuhan-Hu-1"] = []
AA_change_dic["Wuhan-Hu-1"] = {}

"""Load DMS Escape fraction data"""
Escape_Fraction = pd.read_csv(sys.argv[3])
Ab_classes = np.unique(Escape_Fraction["group"].values.astype(str))

"""Load lineage name to assess and it's mutation profile"""
try:
    Lin_name = sys.argv[4]
    mut_x_sites_dic_updated = mut_x_sites_dic.copy()
    AA_change_dic_updated = AA_change_dic.copy()
    if "outbreak" not in Lin_name:
        mut_file = open(sys.argv[5], "r")
        mut_lin0 = mut_file.readlines()
        mut_file.close()
    
        mut_Lin = []
        aa_lin = {}
        for mut in mut_lin0:
            if "\n" in mut:
                mut = mut.replace("\n","")
                
            if mut[:3] not in ("DEL", "del"):
                if len(re.findall(r'\d+', mut))>0:
                    pos0 = re.findall(r'\d+', mut)
                    if len(pos0) == 1:
                        pos = str(pos0[0])
                        mut_Lin.append(pos)   
                        aa_lin[pos] = mut
        """Update mutation profile dictionary"""
        mut_x_sites_dic_updated[Lin_name] = mut_Lin
        AA_change_dic_updated[Lin_name] = aa_lin
except:
    mut_x_sites_dic_updated = mut_x_sites_dic.copy()
    AA_change_dic_updated = AA_change_dic.copy()

"""Relevant functions"""
def get_pos(var_1, var_2, AA_change_dic_1, AA_change_dic_2, mut_x_sites_dic_1, mut_x_sites_dic_2):
    mut_1 = list(AA_change_dic_1[var_1].keys())
    mut_2 = list(AA_change_dic_2[var_2].keys())
    
    pos_diff = list(set(mut_x_sites_dic_1[var_1]).symmetric_difference(set(mut_x_sites_dic_2[var_2]))) 
    sites = []
    if len(mut_1) > len(mut_2):
        for m1 in mut_1:
            for m2 in mut_2:
                if str(m1) == str(m2):
                    check = list(set(AA_change_dic_1[var_1][m1]).intersection(set(AA_change_dic_2[var_2][m2])))
                    if len(check)==0 and (m1 not in sites): # no intersection 
                        sites.append(m1)  
                else:
                    if (m2 not in sites) and (m2 in pos_diff):
                        sites.append(m2)
            
            if (m1 not in sites) and (m1 in pos_diff):
                sites.append(m1)
    else:
        for m2 in mut_2:
            for m1 in mut_1:
                if str(m1) == str(m2):
                    check = list(set(AA_change_dic_1[var_1][m1]).intersection(set(AA_change_dic_2[var_2][m2])))
                    if len(check)==0 and (m2 not in sites): # no intersection 
                        sites.append(m2)  
                else:
                    if (m1 not in sites) and (m1 in pos_diff):
                        sites.append(m1)
            
            if (m2 not in sites) and (m2 in pos_diff):
                sites.append(m2)
    return sites 

def sub_Bind(d, tiled_esc, Where_Mut, Where_Cond):
    Inter_Cond_Mut = Where_Mut & Where_Cond[np.newaxis, d, :]
    Bind_list_d = np.prod(1 - tiled_esc[:, np.newaxis, :]*Inter_Cond_Mut, axis = (1,2))
    Missing_cond_data_d = ~np.any(np.any(Where_Mut & Where_Cond[np.newaxis, d, :], axis = 2), axis = 1)
    return [Bind_list_d, Missing_cond_data_d]


def FR_xy(i, mut_sites, mut_bool_g1, mut_bool_g2, escape_ab_dic, ab, variant_name, aa_diff_bool_i = None, EF_func = "MEAN", GM = False, quiet = True, joblib = None, cluster=False, n_jobs = 10):
    vars_num = mut_bool_g2.shape[0]
    
    test_i = np.tile(mut_bool_g1[i, :], (vars_num, 1))
    
    if aa_diff_bool_i is not None:
        diff_sites = (test_i ^ mut_bool_g2) + aa_diff_bool_i # sites/positions are in symmetric difference and amend FALSE to TRUE if a position is not in the symmetric difference but the aminoacid change is different
    else:
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
    if joblib in (True, "joblib"):
        """Parallel codes --- macOS Monterey 12.5 crashes --- Not used by default """
        pfunc = partial(sub_Bind, tiled_esc = tiled_esc, Where_Mut = Where_Mut, Where_Cond = Where_Cond)
        status = False
        if not cluster:
            try:
                jb_res = list(jb.Parallel(n_jobs = -1, backend = "loky")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
                status = True
                #print("run joblib.Parallel")
            except:
                try:
                    jb_res = list(jb.Parallel(n_jobs = -1, backend = "multiprocessing")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
                    status=True
                    #print("run joblib.Parallel")
                except:
                    jb_res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
                    status=True
                    #print("run joblib.Parallel")
        else:
            #print("run cluster")
            try:
                jb_res = list(jb.Parallel(n_jobs = n_jobs, backend = "multiprocessing")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
                status = True
            except:
                try:
                    jb_res = list(jb.Parallel(n_jobs = n_jobs, backend = "loky")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
                    status=True
                    #print("run joblib.Parallel")
                except:
                    jb_res = list(jb.Parallel(n_jobs = n_jobs, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
                    status=True
        if status:
            for d in range(len(conditions)):
                Bind_list[:, d]   = np.array(jb_res[d][0])
                Missing_cond_data[:, d] = np.array(jb_res[d][1])
        else:
            """ Brute force method """
            print("joblib.Parallel failed running, using brute force looping")
            for d in range(len(conditions)):
                #print(d+1, len(conditions))
                Inter_Cond_Mut = Where_Mut & Where_Cond[np.newaxis, d, :]
                Bind_list[:, d] = np.prod(1 - tiled_esc[:, np.newaxis, :]*Inter_Cond_Mut, axis = (1,2))  
                Missing_cond_data[:, d] = ~np.any(np.any(Where_Mut & Where_Cond[np.newaxis, d, :], axis = 2), axis = 1)
    else:
        """ Brute force method """
        print("Using brute force looping")
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

def cross_reactivity(variant_name, escape_per_sites, Ab_classes, mut_sites_per_variant, AA_change_dic = None, EF_func = "MEAN", GM = False, quiet = True, joblib = None, cluster = False, n_jobs = 10):
    FRxy = {}
    Missed = []
    Greater_one = []
    ### extract variants groups
    variants_g1, variants_g2 = variant_name
    ### extract all sites
    mut_sites = []
    for var in list(mut_sites_per_variant.keys()):
        if (var in variants_g1) or (var in variants_g2):
            mut_sites += list(np.array(mut_sites_per_variant[var]))
    
    mut_sites = np.unique(mut_sites).astype(str)
    ### Construct a boolean array for location of mutation for all variants    
    mut_bool_g1 = np.zeros((len(variants_g1), len(mut_sites)), dtype = bool)
    mut_bool_g2= np.zeros((len(variants_g2), len(mut_sites)), dtype = bool)
    #print(variant_name)
    if AA_change_dic is None:
        aa_bool_diff = None
    else:
        aa_bool_diff = np.zeros((len(variants_g1), len(variants_g2), len(mut_sites))).astype(bool)
        
    for i in range(max(len(variants_g1), len(variants_g2))):
        for j in range(len(mut_sites)): 
            if i<len(variants_g1):
                mut_bool_g1[i, j] = mut_sites[j] in list(np.array(mut_sites_per_variant[variants_g1[i]]).astype(str))
                if AA_change_dic is not None:    
                    for k in range(len(variants_g2)):
                        if (mut_sites[j] in mut_sites_per_variant[variants_g1[i]]) and (mut_sites[j] in mut_sites_per_variant[variants_g2[k]]): ### if the position exists in both variants
                            aa_bool_diff[i, k, j] = len(list(set(AA_change_dic[variants_g1[i]][mut_sites[j]]).intersection(set(AA_change_dic[variants_g2[k]][mut_sites[j]])))) == 0 ### positions will be considered when aa_changes intersection is EMPTY, i.e, all aa changes are different
                
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
            if aa_bool_diff is not None:
                FR, missed, gOne = FR_xy(i, mut_sites, mut_bool_g1, mut_bool_g2, escape_ab_dic, ab, variant_name, aa_diff_bool_i = aa_bool_diff[i, :, :], GM = GM, quiet = quiet, joblib = joblib, cluster = cluster, n_jobs = n_jobs)
            else:
                FR, missed, gOne = FR_xy(i, mut_sites, mut_bool_g1, mut_bool_g2, escape_ab_dic, ab, variant_name, aa_diff_bool_i = None, GM = GM, quiet = quiet, joblib = joblib, cluster = cluster, n_jobs = n_jobs)

            if EF_func == "MEAN":
                FRxy_ab[i, :] = np.mean(FR, axis = 1)
                
            elif EF_func == "PROD":
                FRxy_ab[i, :] = np.prod(FR, axis = 1)
            
            Missing_Cond[i, :] = np.all(missed, axis = 1)
            
            Missed = []
            Greater_one += gOne
            
        #FRxy[ab] = ma.array(FRxy_ab, mask = Missing_Cond.astype(bool), fill_value = np.nan)
        FRxy_ab[Missing_Cond.astype(bool)] = 1
        FRxy[ab] = FRxy_ab.copy()
        
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
        AA_change_dic_updated[spklin] = AA_change_dic[Pseudogroup_dic[spklin]]
    else:           
        try:
            mut_file = open(mut_sim[lineages_sim.index(spklin)], "r")
            mut_lin0 = mut_file.readlines()
            mut_file.close()
            mut_maj= []
            aa_lin = {}
            for mut in mut_lin0:
                if "\n" in mut:
                    mut = mut.replace("\n","")
                    
                if mut[:3] not in ("DEL", "del"):
                    pos0 = re.findall(r'\d+', mut)
                    if len(pos0) == 1:
                        pos = str(pos0[0])
                        if pos not in list(aa_lin.keys()):
                            mut_maj.append(pos)   
                            aa_lin[pos] = [mut]
                        else:
                            aa_lin[pos].append(mut)   
            
            """Update mutation profile dictionary"""        
            mut_x_sites_dic_updated[spklin] = mut_maj
            AA_change_dic_updated[spklin] = aa_lin
            Top_Pseudo.append(spklin)
        except:
            pass

if "Wuhan-Hu-1" not in Top_Pseudo:
    Top_Pseudo = ["Wuhan-Hu-1"] + list(Top_Pseudo)
    
a = 1
if len(Top_Pseudo)!=0:
    Cross_react_dic = {}
    mut_x_sites_dic_used = mut_x_sites_dic_updated.copy()
    AA_change_dic_used = AA_change_dic_updated.copy()
    try:
        mut_x_sites_dic_used[Lin_name] = mut_Lin
        Top_Pseudo.append(Lin_name)
    except:
        pass
      
    # compute mean IC50 per Ab_classes
    IC50_group = Escape_Fraction.groupby('condition', as_index=False).first()[['condition', 'IC50', 'group']]
    mean_IC50_per_group = IC50_group.groupby('group')['IC50'].mean().reset_index()
    print("Mean IC50 per Epitope Classes")
    print(mean_IC50_per_group)
    
    Cross_react_dic_wght = {}
    
    for ab in Ab_classes:
        print("Assess major lineages/pseudogroups with the NTD-RBD mutation positions ")
        print("Cross reactivity countdown", a, "out of %d epitope clases"%len(Ab_classes))
        if ab!= "NTD":
            Cross_Lin, Missed, Greater_one = cross_reactivity((Top_Pseudo, Top_Pseudo), 
                                                                  Escape_Fraction, 
                                                                  [ab],
                                                                  mut_x_sites_dic_used,
                                                                  AA_change_dic = AA_change_dic_updated,
                                                                  joblib=True)
            
            FRxy_ab = Cross_Lin[ab]
            Cross_react_dic[ab] = FRxy_ab.copy()
            Cross_react_dic_wght[ab] = FRxy_ab*((mean_IC50_per_group["IC50"].values[mean_IC50_per_group["group"] == ab])[0])
        a +=1
    
    Cross_react_dic["variant_list"] = list(Top_Pseudo)
    Cross_react_dic_wght["variant_list"] = list(Top_Pseudo)
    
    
    """Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""   
    n = len(Cross_react_dic["variant_list"])
    FR_NTD = np.ones((n, n))
    mut_profiles = []
    for i in range(n):
        var_1 = Cross_react_dic["variant_list"][i]
        if len(list(AA_change_dic_updated[var_1].keys())) > 0:
            var_1_profiles = np.concatenate(tuple([AA_change_dic_updated[var_1][m1] for m1 in list(AA_change_dic_updated[var_1].keys())]))
            mut_profiles.append("/".join(sorted(var_1_profiles)))
        else:
            mut_profiles.append("")
            
        for j in range(n):
            if i > j:
                var_2 = Cross_react_dic["variant_list"][j]
                
                sites = get_pos(var_1, var_2, AA_change_dic_updated, AA_change_dic_updated, mut_x_sites_dic_updated, mut_x_sites_dic_updated)
                
                FR_sites = 1
                pos_done = []
                for s in sites:
                    s = int(s)
                    if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                        if s not in pos_done:
                            FR_sites *= 10
                            pos_done.append(s)

                            
                FR_NTD[i, j] = FR_sites
                FR_NTD[j, i] = FR_sites
        
    Cross_react_dic["NTD"] = FR_NTD.copy()
    Cross_react_dic["Mutations"] = {"mut_profiles":mut_profiles, "positions":mut_x_sites_dic_updated, "AA_changes":AA_change_dic_updated}

    Cross_react_dic_wght["NTD"] = FR_NTD.copy()
    Cross_react_dic_wght["Mutations"] = {"mut_profiles":mut_profiles, "positions":mut_x_sites_dic_updated, "AA_changes":AA_change_dic_updated}

    file0 = open(sys.argv[k], "wb") 
    pickle.dump(Cross_react_dic, file0)
    file0.close()

    file1 = open(str(sys.argv[k])[:-4]+"_weighted.pck", "wb") 
    pickle.dump(Cross_react_dic_wght, file1)
    file1.close()     
