#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 21:19:00 2023

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
try:
    file1 = open(sys.argv[1], "rb") 
    SpikeGroups_list = pickle.load(file1)["names"]
    file1.close()
except:
    SpikeGroups_list = [] ### placeholder for some runs not needing the parameter

"""Load Mutation profile dictionary and aa_changes"""
try:
    file1 = open(sys.argv[2], "rb") 
    Mut_infos = pickle.load(file1)
    mut_x_sites_dic = Mut_infos["positions"]
    AA_change_dic = Mut_infos["AA_changes"]
    file1.close()
except:
    mut_x_sites_dic = {} ### placeholder for some runs not needing the parameter
    AA_change_dic = {} ### placeholder for some runs not needing the parameter

variant_x_names_cross = SpikeGroups_list
# include wild-type
if "Wuhan-Hu-1" not in variant_x_names_cross:
    variant_x_names_cross = ["Wuhan-Hu-1"]+list(variant_x_names_cross)

mut_x_sites_dic["Wuhan-Hu-1"] = []
AA_change_dic["Wuhan-Hu-1"] = {}

"""Load DMS Escape fraction data"""
print(sys.argv[3])
Escape_Fraction = pd.read_csv(sys.argv[3])
Ab_classes = np.unique(Escape_Fraction["group"].values.astype(str))

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
         

"""Compute Cross reactivity"""
Cross_react_dic = {}
AB = ""
a = 1
print("Cross reactivity computation might take a while")
try:
    cluster = str(sys.argv[8])
    if cluster in ("cluster_True", "True", "TRUE"):
        cluster = True
        cluster_argv = True
    else:
        cluster = False
        cluster_argv = False
except:
    cluster = False
    cluster_argv = False
    
try:
    n_jobs = int(sys.argv[9])
    given_njobs = True
except:
    n_jobs = -1
    given_njobs = False

"""Compute Cross reactivity to Delta for valitation"""
if not cluster_argv:
    delta_file = str(sys.argv[-2])
elif cluster_argv and given_njobs:
    delta_file = str(sys.argv[-4])
else:
    delta_file = str(sys.argv[-3])

if delta_file[:4] not in ("None", "none", False):
    # use Delta: B.1.617.2 for validation (Used in Clinical data paper)
    variant_x_names_show = ["Wuhan-Hu-1", "Delta: B.1.617.2"]
    mut_dic_show = {"Wuhan-Hu-1":[], "Delta: B.1.617.2": [614, 950, 142, 452, 681, 19, 478]}
    
    Cross_with_delta_validation, Missed, Greater_one = cross_reactivity((variant_x_names_show,variant_x_names_show), 
    															             Escape_Fraction, Ab_classes, 
                                                            	             mut_dic_show,
                                                                            joblib=True, 
                                                                            cluster = False)
    
    """Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""  
    n = len(variant_x_names_show)
    FR_NTD = np.ones((n, n))
    
    for i in range(n):
        var_1 = variant_x_names_show[i]
        for j in range(n):
            if i > j:
                var_2 = variant_x_names_show[j]
    
                sites_1 = set(np.array(mut_dic_show[var_1]).astype(int))
                sites_2 = set(np.array(mut_dic_show[var_2]).astype(int))
    
                sites = list(sites_1.symmetric_difference(sites_2))
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
    
    Cross_with_delta_validation["NTD"] = FR_NTD.copy()
    Cross_with_delta_validation["variant_list"] = variant_x_names_show
    
    try:
        file0 = open(delta_file, "wb") 
        pickle.dump(Cross_with_delta_validation, file0)
        file0.close()
    except:
        pass

"""Load lineage name to assess and it's mutation profile"""
try:
    n_groups = int(sys.argv[4]) 
    Lin_name = "Groups"
except:
    Lin_name = sys.argv[4]

if Lin_name not in ("ALL", "FR_DMS_sites", "missing", "only_delta"):
    if (Lin_name != "Groups"):
        try:
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
                            if pos not in list(aa_lin.keys()):
                                mut_Lin.append(pos)   
                                aa_lin[pos] = [mut]
                            else:
                                aa_lin[pos].append(mut)

            """Update mutation profile dictionary and AA dictionary """
            mut_x_sites_dic_updated = mut_x_sites_dic.copy()
            AA_change_dic_updated = AA_change_dic.copy()
            if Lin_name not in variant_x_names_cross: ### Keep Lin_name as it is because we advise the user to replace "." with "_" in lineage_focus parameter
                mut_x_sites_dic_updated[Lin_name] = mut_Lin
                AA_change_dic_updated[Lin_name] = aa_lin
            else:
                mut_x_sites_dic_updated[Lin_name] = mut_Lin
                AA_change_dic_updated[Lin_name] = aa_lin
        except:
            sys.exit("Lineage focus mutation file must be provided")
        
        g = []
        g_var =[]
        inds = np.arange(0, len(variant_x_names_cross)).astype(int)
        cut_step = 300
        if len(variant_x_names_cross)>cut_step:
            cut1 = 0
            cut2 = cut_step
            while cut2<len(variant_x_names_cross):
                g.append(inds[cut1:cut2])
                g_var.append(list(np.array(variant_x_names_cross)[cut1:cut2]))
                cut1=cut2
                cut2+=min(cut_step, len(variant_x_names_cross)-cut2)
            g.append(inds[cut1:cut2])
            g_var.append(list(np.array(variant_x_names_cross)[cut1:cut2]))
        else:
            g.append(inds)
            g_var.append(variant_x_names_cross)
        
        if Lin_name not in list(variant_x_names_cross):
            w_lin = len(variant_x_names_cross)
            Cross_react_dic["variant_list"] = list(variant_x_names_cross)+[Lin_name]
        else:
            w_lin = list(variant_x_names_cross).index(Lin_name)
            Cross_react_dic["variant_list"] = list(variant_x_names_cross)
            
        for ab in Ab_classes:
            print("Assess Lineage %s with the NTD-RBD mutation positions "%Lin_name, mut_Lin)
            print("Cross reactivity Epitope %s, countdown"%ab, a, "out of %d epitope clases"%len(Ab_classes))    
            
            if ab!= "NTD":
                
                if Lin_name not in list(variant_x_names_cross):
                    FRxy_ab = np.ones((len(variant_x_names_cross)+1, len(variant_x_names_cross)+1))
                else:
                    FRxy_ab = np.ones((len(variant_x_names_cross), len(variant_x_names_cross)))
                    
                for s in range(len(g)):
                    Cross_Lin, Missed, Greater_one = cross_reactivity(([Lin_name], g_var[s]), 
                                                                              Escape_Fraction, 
                                                                              [ab],
                                                                              mut_x_sites_dic_updated,
                                                                              AA_change_dic=AA_change_dic_updated,
                                                                              joblib=True)
                    
                    #Only the information for the specific lineage studied is required for immunological landscape calculation
                    #the FRxy_ab matrix is kept only for compatibility with other codes
                    
                    FRxy_ab[w_lin, g[s]] = Cross_Lin[ab][0, :]
                    FRxy_ab[g[s], w_lin] = Cross_Lin[ab][0, :]
                Cross_react_dic[ab] = FRxy_ab.copy()
            a +=1
        
        """Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""   
        print("Cross reactivity Epitope NTD")
        n = len(Cross_react_dic["variant_list"])
        FR_NTD = np.ones((n, n))
        mut_profiles = []
        for i1 in range(n):
            var_1 = Cross_react_dic["variant_list"][i1]
            if len(list(AA_change_dic_updated[var_1].keys())) > 0:
                var_1_profiles = np.concatenate(tuple([AA_change_dic_updated[var_1][m1] for m1 in list(AA_change_dic_updated[var_1].keys())]))
                mut_profiles.append("/".join(sorted(var_1_profiles)))
            else:
                mut_profiles.append("")
            for j in range(n):
                if i1 > j:
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
                        
                    FR_NTD[i1, j] = FR_sites
                    FR_NTD[j, i1] = FR_sites
        Cross_react_dic["NTD"] = FR_NTD.copy()
        Cross_react_dic["Mutations"] = {"mut_profiles":mut_profiles, "positions":mut_x_sites_dic_updated, "AA_changes":AA_change_dic_updated}
        file0 = open(sys.argv[7], "wb") 
        pickle.dump(Cross_react_dic, file0)
        file0.close()

    else: ### must be groups
        try:
            file = open("Spikegroups_membership.pck", "rb")
            Pseudogroup_dic = pickle.load(file)
            file.close()
        except:
            Pseudogroup_dic = {} ### placeholder for when parameter for some runs not needing parameter
        
        k = 5
        if (n_groups == 1) and (str(sys.argv[k])[-24:] == "Vaccination_Timeline.csv") and (str(sys.argv[k])[:9] != "outbreak/"): ### Hard coded for vaccination pseudo variants
            vacc_infos = pd.read_csv(sys.argv[k])
            Lin_list = vacc_infos.columns[(vacc_infos.columns != "date")&(vacc_infos.columns != "Unnamed: 0")].tolist()
            mut_sim = ["avail"]*len(Lin_list)
        
        elif (n_groups == 1) and (str(sys.argv[k])[-24:] != "Vaccination_Timeline.csv") and (str(sys.argv[k])[:9] == "outbreak/"): ### Hard coded for oubreak.infos data
            # load special mutation data
            variant_mut_data = pd.read_csv(sys.argv[6])
            Lin_list_0 = str(sys.argv[k])[9:].split("/")
            
            Lin_list = []
            var_outbreak = list(np.array(variant_mut_data["lineage"].values).astype(str))
            mut_outbreak = np.array(variant_mut_data["RBD_NTD_mutations"].values).astype(str)
            
            mut_dic_outbreak = {}
            aa_dic_outbreak = {}
            for lin in Lin_list_0:
                if lin in var_outbreak:
                    w_i = var_outbreak.index(lin)
                    mut_x = mut_outbreak[w_i]
                    split_mut = mut_x.split("/")
                    aa_x = {}
                    pos_list = []
                    for mut in split_mut:
                        pos0 = re.findall(r'\d+', mut)
                        if len(pos0) == 1:
                            pos = str(pos0[0])
                            if pos not in list(aa_x.keys()):
                                aa_x[pos] = [mut]
                                pos_list.append(pos)
                            else:
                                aa_x[pos].append(mut)
                    
                    var_name = "outbreak_%s"%lin
                    mut_dic_outbreak[var_name] = pos_list
                    aa_dic_outbreak[var_name] = aa_x
                    Lin_list.append(var_name)
            
            mut_sim = ["outbreak_placeholder"]*len(Lin_list)
                
        else:
            Lin_list = []
            while k<(n_groups+5):
                Lin_list.append(sys.argv[k])
                k +=1
            
            mut_sim = []
            while k<(n_groups+5+len(Lin_list)):
                mut_sim.append(sys.argv[k])
                k +=1

        mut_x_sites_dic_updated = mut_x_sites_dic.copy()
        AA_change_dic_updated = AA_change_dic.copy()
        Grouped = []
        Lin_exists = []
        Lin_exists_names = []
        single_lin = 0
        Lin_list_grouped = {}
        for j in range(len(Lin_list)):
            Lin_list_grouped[Lin_list[j]] = [Lin_list[j]]
            if Lin_list[j] not in list(Pseudogroup_dic.keys()) and mut_sim[j] != "avail":
                
                if (str(sys.argv[k])[:9] == "outbreak/"): ### Hard coded for oubreak.infos data
                    if j == 0: ## do this only once
                        mut_x_sites_dic_updated.update(mut_dic_outbreak)
                        AA_change_dic_updated.update(aa_dic_outbreak)
                        
                    Grouped.append(False)
                    Lin_exists.append(Lin_list[j])
                    Lin_exists_names.append(Lin_list[j])
                else:
                    try:
                        mut_file = open(mut_sim[j], "r")
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
                                        if pos not in list(aa_lin.keys()):
                                            mut_Lin.append(pos)   
                                            aa_lin[pos] = [mut]
                                        else:
                                            aa_lin[pos].append(mut)   
    
                        """Update mutation profile dictionary"""
                        mut_x_sites_dic_updated[Lin_list[j]] = mut_Lin
                        AA_change_dic_updated[Lin_list[j]] = aa_lin
                        Grouped.append(False)
                        Lin_exists.append(Lin_list[j])
                        Lin_exists_names.append(Lin_list[j])
                    except:
                        try:
                            data_file = open(mut_sim[j], "rb")
                            mutation_loaded = pickle.load(data_file)
                            data_file.close()
                            mutation_data = mutation_loaded["positions"]
                            variants = mutation_loaded["Group"]
                            aa_data = mutation_loaded["AA_changes"]
                            mut_sub = []
                            for i in range(len(variants)):
                                var = variants[i]
                                if var in list(Pseudogroup_dic.keys()):
                                    mut_x_sites_dic_updated[var] = mut_x_sites_dic[Pseudogroup_dic[var]]
                                    AA_change_dic_updated[var] = AA_change_dic[Pseudogroup_dic[var]]
                                else:
                                    mut_x_sites_dic_updated[var] = mutation_data[var]
                                    AA_change_dic_updated[var] = aa_data[var]
                                
                                mut_var = ""
                                for pos in list(AA_change_dic_updated[var].keys):
                                    mut_var += "/".join(AA_change_dic_updated[var][pos])+"/"
                                mut_sub.append(mut_var)
                            
                            inds_spk = []
                            variants_spk = []
                            
                            mut_sub = np.array(mut_sub)
                            for mut in np.unique(mut_sub):
                                var_0 = np.array(variants)[mut_sub == mut][0]
                                variants_spk.append(var_0)
                                inds_spk +=[list(variants).index(var_0)]*np.sum(mut_sub == mut)
    
                            Lin_list_grouped[Lin_list[j]] = variants
                            Lin_list_grouped[Lin_list[j]+"_spk"] = variants_spk
                            Lin_list_grouped["inds"] = inds_spk
                            Grouped.append(True)
                            Lin_exists.append(Lin_list[j])
                            Lin_exists_names.append(Lin_list[j])                       
                            single_lin += len(variants_spk)
                        except:
                            pass
                    
            elif ("*_as_" in Lin_list[j]) and mut_sim[j] == "avail": ### hard-coded in vaccines pseudo names 
                lin_clean = Lin_list[j].split("*_as_")[0]
                mut_x_sites_dic_updated[lin_clean] = mut_x_sites_dic[Pseudogroup_dic[lin_clean]]
                AA_change_dic_updated[lin_clean] = AA_change_dic_updated[Pseudogroup_dic[lin_clean]]
                Grouped.append(False)
                Lin_exists.append(lin_clean)
                Lin_exists_names.append(Lin_list[j])
                
            elif Lin_list[j] in list(Pseudogroup_dic.keys()):
                mut_x_sites_dic_updated[Lin_list[j]] = mut_x_sites_dic[Pseudogroup_dic[Lin_list[j]]]
                AA_change_dic_updated[Lin_list[j]] = AA_change_dic_updated[Pseudogroup_dic[Lin_list[j]]]
                Grouped.append(False)
                Lin_exists.append(Lin_list[j])
                Lin_exists_names.append(Lin_list[j])
                
            else:
                sys.exit("Cannot understand mutation profile informations")
        
        ### Update to available data
        Lin_list = Lin_exists
        Lin_list_names = Lin_exists_names
        
        g = []
        g_var =[]
        inds = np.arange(0, len(variant_x_names_cross)).astype(int)
        
        if single_lin<10:
            cut_step = 300
        else:
            cut_step = 50
        
        if len(variant_x_names_cross)>cut_step:
            cut1 = 0
            cut2 = cut_step
            while cut2<len(variant_x_names_cross):
                g.append(inds[cut1:cut2])
                g_var.append(list(np.array(variant_x_names_cross)[cut1:cut2]))
                cut1=cut2
                cut2+=min(cut_step, len(variant_x_names_cross)-cut2)
            g.append(inds[cut1:cut2])
            g_var.append(list(np.array(variant_x_names_cross)[cut1:cut2]))
        else:
            g.append(inds)
            g_var.append(variant_x_names_cross)
        
        status_sim = []
        for i in range(len(Lin_list)):
            Cross_i = {}
            Cross_i["variant_list"] = list(variant_x_names_cross)
            
            Lin_list_i = Lin_list[i]
            
            if Grouped[i]:
                Lin_list_i = Lin_list_grouped[Lin_list[i]]
                Cross_i["Group"] = Lin_list_i
                
            
            try:
                file_test = open("results/Cross_react_dic_spikegroups_ALL.pck", "rb")
                file_test.close()
                extract = True 
            except:
                extract = False # file is not present and thus if must be recomputed
            
            
            try:
                if (Lin_list[i] not in list(Pseudogroup_dic.keys())) or (not extract):
                    
                    w_lin_i = []
                    add_lin = 0
                    
                    if Lin_list[i]+"_spk" in list(Lin_list_grouped.keys()):
                        Lin_list_i_spk = Lin_list_grouped[Lin_list[i]+"_spk"]
                        inds_spk = np.array(Lin_list_grouped["inds"])
                        spk_adjust = True
                    else:
                        Lin_list_i_spk = Lin_list_grouped[Lin_list[i]]
                        spk_adjust = False
                        inds_spk = np.arange(len(Lin_list_i_spk)).astype(int)
                    
                   
                    for k in range(len(Lin_list_i_spk)):
                        var = Lin_list_i_spk[k]
                        if var not in list(variant_x_names_cross):
                            if not spk_adjust:
                                w_lin_i.append(len(variant_x_names_cross)+add_lin)
                                Cross_i["variant_list"].append(var)
                                add_lin +=1
                            else:
                                where_var = inds_spk == list(Lin_list_i).index(var)
                                w_lin_i += [len(variant_x_names_cross) + j for j in range(add_lin, add_lin+np.sum(where_var))]
                                Cross_i["variant_list"] += list(np.array(Lin_list_i)[where_var])
                                add_lin += np.sum(where_var)
                            
                        else:
                            w_lin = list(variant_x_names_cross).index(var)
                            if not spk_adjust:
                                w_lin_i.append(w_lin)
                            else:
                                where_var = (inds_spk == list(Lin_list_i).index(var))&(np.array(Lin_list_i) != var)
                                w_lin_i += [len(variant_x_names_cross) + j for j in range(add_lin, add_lin+np.sum(where_var))]
                                Cross_i["variant_list"] += list(np.array(Lin_list_i)[where_var])
                                add_lin += np.sum(where_var)
                                        
                    try:
                        file_c = open("results/Cross_react_dic_spikegroups_ALL.pck", "rb") 
                        Cross_global = pickle.load(file_c)
                        variant_global = list(Cross_global["variant_list"])
                        Cross_global.pop("variant_list")
                        file_c.close()
                        Lin_list_i_spk_reduced = []
                        w_global = []
                        w_cross = []
                        for x in Lin_list_i_spk:
                            if x in list(Pseudogroup_dic.keys()):
                                if Pseudogroup_dic[x] in variant_global:
                                    w_global.append(variant_global.index(Pseudogroup_dic[x]))
                                    w_cross.append(list(Cross_i["variant_list"]).index(x))
                                else:
                                    Lin_list_i_spk_reduced.append(x)
                            else:
                                Lin_list_i_spk_reduced.append(x)
                        w_global = np.array(w_global)
                        w_cross = np.array(w_cross)
                    except:
                        Lin_list_i_spk_reduced = Lin_list_i_spk
                        w_global = []
                                        
                    a = 1
                    for ab in Ab_classes:  
                        if ab!= "NTD":
                            if Lin_list[i] not in list(variant_x_names_cross):
                                FRxy_ab = np.ones((len(variant_x_names_cross)+add_lin, len(variant_x_names_cross)+add_lin))
                            else:
                                FRxy_ab = np.ones((len(variant_x_names_cross), len(variant_x_names_cross)))
                                
                            try:
                                print("Assess lineage %s| %d out of %d with the NTD-RBD mutation positions"%(Lin_list[i], i+1,len(Lin_list)), mut_x_sites_dic_updated[Lin_list[i]])
                            except:
                                print("Assess lineage %s (%d spikesgroups, %d lineages)| %d out of %d with the NTD-RBD mutation positions"%(Lin_list[i], len(Lin_list_i_spk), len(Lin_list_i), i+1,len(Lin_list)))
                            
                            print("Cross reactivity Epitope %s, countdown"%ab, a, "out of %d epitope clases"%len(Ab_classes)) 
                            sub_FR = np.ones((len(Cross_i["variant_list"]), len(Cross_i["variant_list"])))
                            
                            if len(w_global)>0:
                                for ex in range(len(w_cross)):
                                    sub_FR[w_cross[ex], w_cross] = Cross_global[ab][w_global[ex], w_global]
                            
                                    
                            if len(Lin_list_i_spk_reduced)>0:
                                where_spk_s = np.array([list(Cross_i["variant_list"]).index(Lin_list_i_spk_reduced[k]) for k in range(len(Lin_list_i_spk_reduced))]) 
                                for s in range(len(g)):
                                    Cross_Lin, Missed, Greater_one = cross_reactivity((Lin_list_i_spk_reduced, g_var[s]), 
                                                                                      Escape_Fraction, 
                                                                                      [ab],
                                                                                      mut_x_sites_dic_updated,
                                                                                      AA_change_dic=AA_change_dic_updated,
                                                                                      joblib=True)
                                
                               
                                    #Only the information for the specific lineage studied is required for immunological landscape calculation
                                    #the FRxy_ab matrix is kept only for compatibility with other codes
                                    locs = np.array([list(Cross_i["variant_list"]).index(g_var[s][k]) for k in range(len(g[s]))])
                                    
                                    for k in range(len(Lin_list_i_spk_reduced)):
                                        sub_FR[where_spk_s[k], locs] = Cross_Lin[ab][k, :]
                                        sub_FR[locs, where_spk_s[k]] = Cross_Lin[ab][k, :]
                            
                                if len(Lin_list_i_spk_reduced)>1:
                                    for k in range(len(Lin_list_i_spk_reduced)):
                                        Cross_Lin, Missed, Greater_one = cross_reactivity(([Lin_list_i_spk_reduced[k]], Lin_list_i_spk_reduced[k+1:]), 
                                                                                  Escape_Fraction, 
                                                                                  [ab],
                                                                                  mut_x_sites_dic_updated,
                                                                                  AA_change_dic=AA_change_dic_updated,
                                                                                  joblib=True)
                                
                                        sub_FR[where_spk_s[k], where_spk_s[k+1:]] = Cross_Lin[ab][k, :]
                                        sub_FR[where_spk_s[k+1:], where_spk_s[k]] = sub_FR[where_spk_s[k], where_spk_s[k+1:]]
                            
                            if not spk_adjust:
                                FRxy_ab[np.array(w_lin_i), :] = sub_FR[np.array(w_lin_i), :]
                                if len(w_lin_i)==1:
                                    FRxy_ab[:, np.array(w_lin_i)] = sub_FR[:, np.array(w_lin_i)]
                                else:
                                    FRxy_ab[:, np.array(w_lin_i)] = sub_FR[:, np.array(w_lin_i)].T
                            else:
                                for w_spk in range(len(Lin_list_i_spk)):
                                    w_spk_cross = list(Cross_i["variant_list"]).index(Lin_list_i_spk[w_spk])
                                    where_spk = inds_spk == list(Lin_list_i).index(Lin_list_i_spk[w_spk])
                                    
                                    Lins = np.array(Lin_list_i)[where_spk]
                                    cross_spk = np.row_stack(tuple([sub_FR[w_spk_cross, :]]*np.sum(where_spk)))
                                    
                                    where_spk_cross = np.array([list(Cross_i["variant_list"]).index(Lins[k]) for k in range(len(Lins))]) 
                                    FRxy_ab[where_spk_cross, :] = cross_spk
                                    if len(where_spk_cross)>1:
                                        FRxy_ab[:, where_spk_cross] = cross_spk.T
                                    else:
                                        FRxy_ab[:, where_spk_cross] = cross_spk
                                        
                            Cross_i[ab] = FRxy_ab.copy()
                        a +=1 
                    
                    """Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""
                    print("Cross reactivity Epitope NTD")
                    n = len(Cross_i["variant_list"])
                    FR_NTD = np.ones((n, n))
                    for i1 in range(n):
                        var_1 = Cross_i["variant_list"][i1]
                        for j1 in range(n):
                            if i1 > j1:
                                var_2 = Cross_i["variant_list"][j1]
                    
                                sites = get_pos(var_1, var_2, AA_change_dic_updated, AA_change_dic_updated, mut_x_sites_dic_updated, mut_x_sites_dic_updated)
                                
                                FR_sites = 1
                                pos_done = []
                                for s in sites:
                                    s = int(s)
                                    if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                                        if s not in pos_done:
                                            FR_sites *= 10
                                            pos_done.append(s)

                                FR_NTD[i1, j1] = FR_sites
                                FR_NTD[j1, i1] = FR_sites
                    Cross_i["NTD"] = FR_NTD.copy()
    
                else:
                    ### open global cross_reactivity file which must be present
                    a = 1  
                    print("Assess lineage %s| %d out of %d with the NTD-RBD mutation positions"%(Lin_list[i], i+1,len(Lin_list)), mut_x_sites_dic_updated[Lin_list[i]])
                    print("Load : %s is present in general file results/Cross_react_dic_spikegroups_ALL.pck"%Lin_list[i]) 
                
                    file_c = open("results/Cross_react_dic_spikegroups_ALL.pck", "rb") 
                    Cross_global = pickle.load(file_c)
                    variant_global = Cross_global["variant_list"]
                    Cross_global.pop("variant_list")
                    try:
                        Cross_global.pop("Mutations")
                    except:
                        pass

                    Ab_global = Cross_global.keys()
                    file_c.close()

                    w_lin = len(variant_x_names_cross)
                    
                    Cross_i["variant_list"] = list(variant_x_names_cross)
                    ### Avoid that Lin_list[i] is sited multiple times Cross_i["variant_list"] because the other entries are not computed
                    if Lin_list[i] in list(Pseudogroup_dic.keys()):
                        if Lin_list[i] == Pseudogroup_dic[Lin_list[i]]:
                        	Cross_i["variant_list"][list(Cross_i["variant_list"]).index(Lin_list[i])] = Lin_list[i]+"_ignore"
                    
                    Cross_i["variant_list"] = Cross_i["variant_list"] + [Lin_list[i]]
                    for ab in Ab_global:
                        FRxy_ab = np.ones((len(variant_x_names_cross)+1, len(variant_x_names_cross)+1))
                            
                        for u1 in range(len(variant_x_names_cross)):
                            v_u1 = variant_x_names_cross[u1]
                            if v_u1 in variant_global:
                                FRxy_ab[w_lin, u1] = Cross_global[ab][list(variant_global).index(Pseudogroup_dic[Lin_list[i]]), list(variant_global).index(v_u1)]
                                FRxy_ab[u1, w_lin] = FRxy_ab[w_lin, u1]
                            else:
                                # recompute it, should not happen normaly
                                # usefull when using previously computed cross reactivity file where all covsonar lineages (not only spikegroups) are present (assigned the FR of their spikegroups)
                                if ab != "NTD":
                                    Cross_Lin, Missed, Greater_one = cross_reactivity(([Lin_list[i]], [v_u1]), 
                                                                                          Escape_Fraction, 
                                                                                          [ab],
                                                                                          mut_x_sites_dic_updated,
                                                                                          AA_change_dic= AA_change_dic_updated,
                                                                                          joblib=True)
                                    
                                   
                                    #Only the information for the specific lineage studied is required for immunological landscape calculation
                                    #the FRxy_ab matrix is kept only for compatibility with other codes
                                    FRxy_ab[w_lin, u1] = Cross_Lin[ab][0, 0]
                                    FRxy_ab[u1, w_lin] = FRxy_ab[w_lin, u1]
                                else:
                                    var_1 = Lin_list[i]
                                    var_2 = v_u1
        
                                    sites = get_pos(var_1, var_2, AA_change_dic_updated, AA_change_dic_updated, mut_x_sites_dic_updated, mut_x_sites_dic_updated)
                                    
                                    FR_sites = 1
                                    pos_done = []
                                    for s in sites:
                                        s = int(s)
                                        if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                                            if s not in pos_done:
                                                FR_sites *= 10
                                                pos_done.append(s)
                                    
                                                    
                                    FRxy_ab[w_lin, u1] = FR_sites
                                    FRxy_ab[u1, w_lin] = FR_sites
                                                       
                        
                        Cross_i[ab] = FRxy_ab.copy()
                        a +=1 
                
                status_sim.append("Done")
                
                n = len(Cross_i["variant_list"])

                mut_profiles = []
                for i1 in range(n):
                    var_1 = Cross_i["variant_list"][i1]
                    if var_1[-len("_ignore"):] == "_ignore":
                        var_1a = var_1[:-len("_ignore")]
                        if len(list(AA_change_dic_updated[var_1a].keys())) > 0:
                            var_1_profiles = np.concatenate(tuple([AA_change_dic_updated[var_1a][m1] for m1 in list(AA_change_dic_updated[var_1a].keys())]))
                            mut_profiles.append("/".join(sorted(var_1_profiles)))
                        else:
                            mut_profiles.append("")
                    else:
                        if len(list(AA_change_dic_updated[var_1].keys())) > 0:
                            var_1_profiles = np.concatenate(tuple([AA_change_dic_updated[var_1][m1] for m1 in list(AA_change_dic_updated[var_1].keys())]))
                            mut_profiles.append("/".join(sorted(var_1_profiles)))
                        else:
                            mut_profiles.append("")
                
                          
                Cross_i["Mutations"] = {"mut_profiles":mut_profiles, "positions":mut_x_sites_dic_updated, "AA_changes":AA_change_dic_updated}
                file0 = open(sys.argv[len(sys.argv)-1]+"/Cross_%s.pck"%Lin_list_names[i], "wb") 
                pickle.dump(Cross_i, file0)
                file0.close()
            except:
                status_sim.append("Failed")
                print("Ignored Cross of %s: Give mutation file or chose a lineage with pseudogroup present in covsonar data"%Lin_list_names[i])
                status_sim.append("No mutation profile for %s, give mutation file or chose a lineage present in covsonar data"%Lin_list_names[i])
            
        if len(status_sim)>0:
            stat_df = pd.DataFrame({"Lineages":Lin_list_names, "computed_cross":status_sim})
            stat_df.to_csv(sys.argv[len(sys.argv)-1]+"/cross_status.csv")
            
elif Lin_name == "ALL":  
    """Break runs into manageable pieces"""
    
    g = []
    g_var =[]
    inds = np.arange(0, len(variant_x_names_cross)).astype(int)
    
    if len(variant_x_names_cross)>100:
        print("WARNING: Cross reactivity computation might take really long on local computers, better using HPC (see scripts/run_cross.sh for hints on using a slurm cluster)")
    
    cut_step = 200
    if len(variant_x_names_cross)>cut_step:
        cut1 = 0
        cut2 = cut_step
        while cut2<len(variant_x_names_cross):
            g.append(inds[cut1:cut2])
            g_var.append(list(np.array(variant_x_names_cross)[cut1:cut2]))
            cut1=cut2
            cut2+=min(cut_step, len(variant_x_names_cross)-cut2)
        g.append(inds[cut1:cut2])
        g_var.append(list(np.array(variant_x_names_cross)[cut1:cut2]))
    else:
        g.append(inds)
        g_var.append(variant_x_names_cross)
    
    # FR is symmetric, computing the upper diagonal or the lower diagonal is enough, thus reducing compuation time to half
    for ab in Ab_classes:
        if ab!= "NTD":
            FRxy_ab = np.zeros((len(variant_x_names_cross), len(variant_x_names_cross)))
            for s1 in range(len(g)):
                print("Assess all spikegroups with the NTD-RBD mutation positions ")
                print("Cross reactivity Epitope %s, countdown"%ab, a, "out of %d epitope clases"%len(Ab_classes))
                print("Run Cross for %d out of %d (to achieve %d/%d spikegroups)"%(s1+1, len(g), len(g[s1]), len(variant_x_names_cross)))
                sub_FR = np.zeros((len(g[s1]), len(variant_x_names_cross)))
                
                ## deal with diagonal block
                for j in range(len(g[s1])):
                    var_sub = np.array(variant_x_names_cross)[g[s1][j]+1:]
                    inds_sub = inds[g[s1][j]+1:]
                    g2 = []
                    g2_var = []
                    if len(var_sub)>cut_step:
                        cut1 = 0
                        cut2 = cut_step
                        while cut2<len(var_sub):
                            g2.append(inds_sub[cut1:cut2])
                            g2_var.append(list(var_sub[cut1:cut2]))
                            cut1 = cut2
                            cut2 += min(cut_step, len(var_sub) - cut2)
                        g2.append(inds_sub[cut1:cut2])
                        g2_var.append(list(var_sub[cut1:cut2]))
                    else:
                        g2.append(inds_sub)
                        g2_var.append(var_sub)
                    
                    for s2 in range(len(g2)):
                        Cross_Lin, Missed, Greater_one = cross_reactivity(([g_var[s1][j]], g2_var[s2]), 
                                                                       Escape_Fraction, 
                                                                       [ab],
                                                                       mut_x_sites_dic, 
                                                                       AA_change_dic = AA_change_dic,
                                                                       joblib=True, 
                                                                       cluster = cluster,
                                                                       n_jobs = n_jobs)
    
                        sub_FR[j, g2[s2]] = Cross_Lin[ab][0, :]
                    
                # deal with off-diagonal block
                var_sub = np.array(variant_x_names_cross)[g[s1][-1]+1:]
                inds_sub = inds[g[s1][-1]+1:]
                g2 = []
                g2_var = []
                if len(var_sub)>cut_step:
                    cut1 = 0
                    cut2 = cut_step
                    while cut2<len(var_sub):
                        g2.append(inds_sub[cut1:cut2])
                        g2_var.append(list(var_sub[cut1:cut2]))
                        cut1 = cut2
                        cut2 += min(cut_step, len(var_sub) - cut2)
                    g2.append(inds_sub[cut1:cut2])
                    g2_var.append(list(var_sub[cut1:cut2]))
                else:
                    g2.append(inds_sub)
                    g2_var.append(var_sub)
                
                for s2 in range(len(g2)):
                    Cross_Lin, Missed, Greater_one = cross_reactivity((g_var[s1], g2_var[s2]), 
                                                                   Escape_Fraction, 
                                                                   [ab],
                                                                   mut_x_sites_dic, 
                                                                   AA_change_dic = AA_change_dic,
                                                                   joblib=True, 
                                                                   cluster = cluster,
                                                                   n_jobs = n_jobs)

                    sub_FR[:, g2[s2]] = Cross_Lin[ab]
                
                FRxy_ab[g[s1], :] = sub_FR
            
            Cross_react_dic[ab] = FRxy_ab + FRxy_ab.T + np.diag(np.ones(len(variant_x_names_cross)))
        a +=1
        
    Cross_react_dic["variant_list"] = list(variant_x_names_cross)
    """Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""  
    print("Cross reactivity spikegroups for Epitope NTD")
    n = len(Cross_react_dic["variant_list"])
    FR_NTD = np.ones((n, n))
    mut_profiles = []
    for i in range(n):
        var_1 = Cross_react_dic["variant_list"][i]
        if len(list(AA_change_dic[var_1].keys())) > 0:
            var_1_profiles = np.concatenate(tuple([AA_change_dic[var_1][m1] for m1 in list(AA_change_dic[var_1].keys())]))
            mut_profiles.append("/".join(sorted(var_1_profiles)))
        else:
            mut_profiles.append("")
        var_1 = Cross_react_dic["variant_list"][i]
        for j in range(n):
            if i > j:
                var_2 = Cross_react_dic["variant_list"][j]
    
                sites = get_pos(var_1, var_2, AA_change_dic, AA_change_dic, mut_x_sites_dic, mut_x_sites_dic)
                
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
    Cross_react_dic["Mutations"] = {"mut_profiles":mut_profiles, "positions":mut_x_sites_dic, "AA_changes":AA_change_dic}
    file0 = open(sys.argv[7], "wb") 
    pickle.dump(Cross_react_dic, file0)
    file0.close()

elif Lin_name == "missing":
    file_c = open(sys.argv[5], "rb") ### hard-coded specifically for cross_missing = TRUE, sys.argv[5] = results/Cross_react_dic_spikegroups_ALL.pck
    Cross_global = pickle.load(file_c)
    variant_global = list(Cross_global["variant_list"])
    
    Cross_global.pop("variant_list")
    
    mut_x_global = Cross_global["Mutations"]["positions"]
    AA_global = Cross_global["Mutations"]["AA_changes"]
    mut_profiles_global = Cross_global["Mutations"]["mut_profiles"]
    
    Cross_global.pop("Mutations")
    
    Ab_global = Cross_global.keys()
    file_c.close()
    """Mutation profiles of variant_x_names_cross"""
    profiles_x_names_cross = []
    for vn in range(len(variant_x_names_cross)):
        lin  = variant_x_names_cross[vn]
        if len(list(AA_change_dic[lin].keys())) > 0:
            var_1_profiles = np.concatenate(tuple([AA_change_dic[lin][m1] for m1 in list(AA_change_dic[lin].keys())]))
            lin_profile = "/".join(sorted(var_1_profiles))
        else:
            lin_profile = ""
        profiles_x_names_cross.append(lin_profile)
    
    """Find indexes of missing and not missing variants"""
    Lin_miss = []
    loc_not_miss = []
    loc_in_cross = []
    sub_miss = {}
    num_rerun = []
    for vn in range(len(variant_x_names_cross)):
        lin  = variant_x_names_cross[vn]
        lin_profile = profiles_x_names_cross[vn]
        if lin_profile not in mut_profiles_global:
            Lin_miss.append(lin)
            sub_miss[lin] = np.ones(len(variant_global)).astype(bool)
            num_rerun.append(len(variant_global))
        else:
            indx = list(mut_profiles_global).index(lin_profile)
            var_1 = variant_global[indx]
            var_2 = lin
            
            #if var_1 != var_2: 
            #    print("diff lineages", var_1, var_2, lin_profile)
                
            sites = get_pos(var_1, var_2, AA_global, AA_change_dic, mut_x_global, mut_x_sites_dic)
            
            if (len(sites) != 0):
                """This part still needs to be proof-checked: normaly this never occurs since lin_profile is the same for var_1 and var_2 (should be refined and added to the case when lin_profile not in mut_profile_global)"""
                # mutation profile is different from general file, thus must be recomputed
                sub_miss[lin] = np.ones(len(variant_global)).astype(bool)
                for ig in range(len(variant_global)):
                    if len(list(AA_global[variant_global[ig]].keys())) > 0:
                        var_1_profiles = np.concatenate(tuple([AA_global[variant_global[ig]][m1] for m1 in list(AA_global[variant_global[ig]].keys())]))
                        ig_profile = "/".join(sorted(var_1_profiles))
                    else:
                        ig_profile = ""

                    if ig_profile in profiles_x_names_cross:
                        ig_in_cross = variant_x_names_cross[profiles_x_names_cross.index(ig_profile)]
                        sites_x = get_pos(ig_in_cross, lin, AA_change_dic, AA_change_dic, mut_x_sites_dic, mut_x_sites_dic)
                        sites_g = get_pos(variant_global[ig], var_1, AA_global, AA_global, mut_x_global, mut_x_global)
                        "check if the symmetric difference in Cross_global is the same as the symmetric difference that we need to consider"
                        if "/".join(sorted(sites_x)) == "/".join(sorted(sites_g)): 
                            sub_miss[lin][ig] = False
                
                Lin_miss.append(lin)
            else:
                loc_in_cross.append(list(variant_x_names_cross).index(var_2))
                loc_not_miss.append(list(variant_global).index(var_1))
                sub_miss[lin] = np.zeros(len(variant_global)).astype(bool)
                
            num_rerun.append(np.sum(sub_miss[lin]))
    
    if len(Lin_miss) == 0:
        Cross_react_dic = Cross_global.copy()
        Cross_react_dic["variant_list"] = list(np.array(variant_x_names_cross)[np.array(loc_in_cross)])
        n = len(Cross_react_dic["variant_list"])
        FR_NTD = np.ones((n, n))
        mut_profiles = []
        for i in range(n):
            var_1 = Cross_react_dic["variant_list"][i]
            if len(list(AA_change_dic[var_1].keys())) > 0:
                var_1_profiles = np.concatenate(tuple([AA_change_dic[var_1][m1] for m1 in list(AA_change_dic[var_1].keys())]))
                mut_profiles.append("/".join(sorted(var_1_profiles)))
            else:
                mut_profiles.append("")
            
            for j in range(n):
                if i > j:
                    var_2 = Cross_react_dic["variant_list"][j]
                    
                    sites = get_pos(var_1, var_2, AA_change_dic, AA_change_dic, mut_x_sites_dic, mut_x_sites_dic)
                    
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
    else:     
        av_rerun = np.mean(num_rerun)
        if len(Lin_miss)>10 and av_rerun > 100:
            print("WARNING: Cross reactivity computation might take really long on local computers, better using HPC (see scripts/run_cross_missing.sh for hints on using a slurm cluster)")
        
        if len(loc_in_cross)!=0:
            w_in_cross = np.arange(0, len(variant_x_names_cross)).astype(int)[np.array(loc_in_cross)]
            variants_in_global = np.array(variant_x_names_cross)[w_in_cross]
            Cross_react_dic["variant_list"] = list(variants_in_global) + Lin_miss
            variants_to_sim = variants_in_global.copy()
        else:
            w_in_cross = np.arange(0, len(variant_x_names_cross)).astype(int)
            variants_to_sim = np.array(variant_x_names_cross)[w_in_cross]
            Cross_react_dic["variant_list"] = list(variants_to_sim)
            
        a = 1
        g = []
        g_var =[]
        inds = np.arange(0, len(variants_to_sim)).astype(int)
        cut_step = 100
        if len(variants_to_sim)>cut_step:
            cut1 = 0
            cut2 = cut_step
            while cut2<len(variants_to_sim):
                g.append(inds[cut1:cut2])
                g_var.append(list(np.array(variants_to_sim)[cut1:cut2]))
                cut1=cut2
                cut2+=min(cut_step, len(variants_to_sim)-cut2)
            g.append(inds[cut1:cut2])
            g_var.append(list(np.array(variants_to_sim)[cut1:cut2]))
        else:
            g.append(inds)
            g_var.append(variants_to_sim)
        
        if len(loc_not_miss) != 0:
            w_global = np.arange(0, len(variant_global)).astype(int)[np.array(loc_not_miss)] ## location of variant_x_names cross in global file
        else:
            w_global = None
        
        # get mut profiles of missing lineages
        lin_profile_list = []
        for pr in range(len(Lin_miss)):
            if len(list(AA_change_dic[lin].keys())) > 0:
                var_1_profiles = np.concatenate(tuple([AA_change_dic[lin][m1] for m1 in list(AA_change_dic[lin].keys())]))
                lin_profile = "/".join(sorted(var_1_profiles))
            else:
                lin_profile = ""
            lin_profile_list.append(lin_profile)
        
        mut_profiles = []
        for i in range(len(Cross_react_dic["variant_list"])):
            var_1 = Cross_react_dic["variant_list"][i]
            if len(list(AA_change_dic[var_1].keys())) > 0:
                var_1_profiles = np.concatenate(tuple([AA_change_dic[var_1][m1] for m1 in list(AA_change_dic[var_1].keys())]))
                mut_profiles.append("/".join(sorted(var_1_profiles)))
            else:
                mut_profiles.append("")
        
        for ab in Ab_global:
            Cross_react_dic[ab] = np.ones((len(Cross_react_dic["variant_list"]), len(Cross_react_dic["variant_list"])))
            if (ab != "NTD") and (ab in Ab_classes):
                for indx_lin in range(len(Lin_miss)):
                    lin = Lin_miss[indx_lin]
                    lin_indx = list(Cross_react_dic["variant_list"]).index(lin)
                    
                    lin_profile = lin_profile_list[indx_lin]
                    
                    n_0 = 0
                    for s in range(len(g)):
                        recomp_lin = [k for k in range(len(g[s])) if sub_miss[lin][mut_profiles_global.index(profiles_x_names_cross[variant_x_names_cross.index(g_var[s][k])])] if profiles_x_names_cross[variant_x_names_cross.index(g_var[s][k])] in mut_profiles_global]
                        recomp_lin += [k for k in range(len(g[s])) if profiles_x_names_cross[variant_x_names_cross.index(g_var[s][k])] not in mut_profiles_global]

                        g_var_recompute = np.array(g_var[s])[np.array(recomp_lin)]
                        
                        if len(g_var_recompute) != 0:
                            print("Assess missing | num %d-th out of %d, to recompute (%d, %d (max %d)) with the NTD-RBD mutation positions"%(indx_lin+1, len(Lin_miss), len(g_var_recompute), min(n_0 + len(g_var_recompute), len(variants_to_sim)),len(variants_to_sim)))
                            print("Cross reactivity Epitope %s, countdown"%ab, a, "out of %d epitope clases"%len(Ab_global))
                            Cross_Lin, Missed, Greater_one = cross_reactivity(([lin], g_var_recompute), 
                                                                              Escape_Fraction, 
                                                                              [ab],
                                                                              mut_x_sites_dic,
                                                                              AA_change_dic = AA_change_dic,
                                                                              joblib = True,
                                                                              cluster = cluster,
                                                                              n_jobs = n_jobs)
                            
                            locs_recompute = np.array([list(Cross_react_dic["variant_list"]).index(g_var_recompute[i]) for i in range(len(g_var_recompute))])
                            Cross_react_dic[ab][lin_indx, locs_recompute] = Cross_Lin[ab]
                            Cross_react_dic[ab][locs_recompute, lin_indx] = Cross_Lin[ab]
                            
                        if len(g_var_recompute) != len(g_var[s]):
                            """This part still needs to be proof-checked: not occuring at the moment (needs to be refined) """
                            if lin_profile in mut_profiles_global:
                                id_lin_global = list(mut_profiles_global).index(lin_profile)
                                g_not_recomputed = [g_var[s][k] for k in range(len(g_var[s])) if g_var[s][k] not in list(g_var_recompute)]                                    
                                if len(g_not_recomputed) != 0:
                                    locs_not_recompt = np.array([list(Cross_react_dic["variant_list"]).index(g_not_recomputed[i]) for i in range(len(g_not_recomputed))])
                                    keep = np.array([mut_profiles_global.index(profiles_x_names_cross[variant_x_names_cross.index(g_not_recomputed[i])]) for i in range(len(g_not_recomputed)) if profiles_x_names_cross[variant_x_names_cross.index(g_not_recomputed[i])] in mut_profiles_global])
                                    if len(keep)>0:
                                        Cross_react_dic[ab][lin_indx, locs_not_recompt] = Cross_global[ab][id_lin_global, keep]
                                        Cross_react_dic[ab][locs_not_recompt, lin_indx] = Cross_global[ab][keep, id_lin_global]
                        
                        n_0 += len(g_var_recompute)
                            
                    print("Assess %d-th missing vs. %d missing with the NTD-RBD mutation positions"%(indx_lin+1, len(Lin_miss)))
                    sub_miss_reduced = [i for i in range(len(Lin_miss)) if (lin_profile_list[i] in mut_profiles_global) and (sub_miss[lin][list(mut_profiles_global).index(lin_profile_list[i])])]
                    sub_miss_reduced += [i for i in range(len(Lin_miss)) if (lin_profile_list[i] not in mut_profiles_global)]
                    recomp_lin_miss = np.array(sub_miss_reduced)
                    Lin_miss_recompute = list(np.array(Lin_miss)[recomp_lin_miss])
                    
                    if len(Lin_miss_recompute)>0:
                        Cross_Lin, Missed, Greater_one = cross_reactivity(([lin], Lin_miss_recompute), 
                                                                  Escape_Fraction, 
                                                                  [ab],
                                                                  mut_x_sites_dic,
                                                                  AA_change_dic = AA_change_dic,
                                                                  joblib = True,
                                                                  cluster = cluster,
                                                                  n_jobs = n_jobs)
                        
                        miss_locs = np.array([list(Cross_react_dic["variant_list"]).index(Lin_miss_recompute[i]) for i in range(len(Lin_miss_recompute))])
                        Cross_react_dic[ab][lin_indx, miss_locs] = Cross_Lin[ab]
                        Cross_react_dic[ab][miss_locs, lin_indx] = Cross_Lin[ab]
                    
                    if len(Lin_miss_recompute)!= len(Lin_miss):
                        if lin_profile in mut_profiles_global:
                            """This part still needs to be proof-checked: not occuring at the moment (needs to be refined) """
                            id_lin_global = list(mut_profiles_global).index(lin_profile)
                            Lin_miss_not_recomputed = []
                            not_recomputed = [i for i in range(len(Lin_miss)) if (lin_profile_list[i] in mut_profiles_global) and not (sub_miss[lin][list(mut_profiles_global).index(lin_profile_list[i])])]
                            present_indx = np.array(not_recomputed)
                            if len(present_indx)>0:
                                Lin_miss_not_recomputed = list(np.array(Lin_miss)[present_indx])
                                profile_not_recomputed = list(np.arry(lin_profile_list)[present_indx])
                                for ind2 in range(len(Lin_miss_not_recomputed)):
                                    loc2 = list(Cross_react_dic["variant_list"]).index(Lin_miss_not_recomputed[ind2])
                                    ind2_global = list(mut_profiles_global).index(profile_not_recomputed[ind2])
                                    Cross_react_dic[ab][lin_indx, loc2] = Cross_global[ab][id_lin_global, ind2_global]
                                    Cross_react_dic[ab][loc2, lin_indx] = Cross_global[ab][ind2_global, id_lin_global]
                
                if w_global is not None:
                    ### Include not missing and not recomputed
                    for lin in Cross_react_dic["variant_list"]:
                        id_lin = list(Cross_react_dic["variant_list"]).index(lin)
                        lin_profile = mut_profiles[id_lin]
                        if lin_profile in mut_profiles_global:
                            id_lin_global = list(mut_profiles_global).index(lin_profile)
                            not_miss_recomputed = w_global[~sub_miss[lin][np.array(loc_not_miss)]]
                            Cross_react_dic[ab][:, :len(variants_in_global)][id_lin, :] = Cross_global[ab][id_lin_global, not_miss_recomputed]
                            Cross_react_dic[ab][:len(variants_in_global), :][:, id_lin] = Cross_global[ab][not_miss_recomputed, id_lin_global]
                            
                a +=1
        
        print("Cross reactivity spikegroups for Epitope NTD")
        n = len(Cross_react_dic["variant_list"])
        FR_NTD = np.ones((n, n))
        for i in range(n):
            var_1 = Cross_react_dic["variant_list"][i]
            for j in range(n):
                if i > j:
                    var_2 = Cross_react_dic["variant_list"][j]
                    
                    sites = get_pos(var_1, var_2, AA_change_dic, AA_change_dic, mut_x_sites_dic, mut_x_sites_dic)
                    
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
            
    Cross_react_dic["Mutations"] = {"mut_profiles":mut_profiles, "positions":mut_x_sites_dic, "AA_changes":AA_change_dic}
    file0 = open(sys.argv[7], "wb") 
    pickle.dump(Cross_react_dic, file0)
    file0.close()

elif Lin_name == "FR_DMS_sites":
    """ Compute FR sites DMS """
    One_mut_lin = np.unique(Escape_Fraction["site"].values.astype(str))
    """Add NTD sites """
    
    One_mut_lin = np.concatenate((One_mut_lin, np.arange(14, 21), np.arange(140, 159), np.arange(245, 265))).astype(int)
    One_mut_lin = np.unique(np.array(One_mut_lin))
    One_mut_lin = One_mut_lin[np.argsort(One_mut_lin)]
    One_mut_lin = One_mut_lin.astype(str)
    
    Ab_One_Mut = Ab_classes
    One_mut_dic = {}
    for x in One_mut_lin:
        One_mut_dic[x] = [x]
    
    """ Add a wild-type lineage place holder"""
    One_mut_lin_new = np.append(["WT"], One_mut_lin)
    One_mut_dic["WT"] = []
    
    """Compute FR for each sites and Ab"""
    FR_Sites_Ab = np.ones((len(Ab_One_Mut), len(One_mut_lin_new)))
    for k in range(len(Ab_One_Mut)):
        ab = Ab_One_Mut[k]
        print("Cross reactivity DMS sites countdown %d out of %d epitope clases"%(k+1, len(Ab_One_Mut)))
        FR, missed, gOne = cross_reactivity((["WT"], One_mut_lin_new),
                                                Escape_Fraction, 
                                                [ab],
                                                One_mut_dic,
                                                joblib=True,
                                                cluster = False)
        
        FR_Sites_Ab[k, :] = FR[ab][0, :]
      
    """Add NTD-targeting antibody class"""
    Ab_One_Mut = list(Ab_One_Mut) + ["NTD"]
    
    """
    Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite
    Reference to Antigenic Supersites: McCallum et al. 2021: N-terminal domain antigenic mapping reveals a site 
    vulnerability of SARS-Cov-2
    """   
    idx_WT = list(One_mut_lin_new).index("WT")
    n = len(One_mut_lin_new)
    FR_NTD = np.ones(n)
    for i in range(n):
        var_1 = One_mut_lin_new[i]
        if i > idx_WT:
            var_2 = One_mut_lin_new[idx_WT]
    
            sites_1 = set(np.array(One_mut_dic[var_1]).astype(int))
            sites_2 = set(np.array(One_mut_dic[var_2]).astype(int))
    
            sites = list(sites_1.symmetric_difference(sites_2))
            FR_sites = 1
            pos_done = []
            for s in sites:
                s = int(s)
                if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                    if s not in pos_done:
                        FR_sites *= 10
                        pos_done.append(s)
            FR_NTD[i] = FR_sites
    
    FR_Sites_Ab = np.row_stack((FR_Sites_Ab, FR_NTD)) 
    
    
    ### Saving file
    FR_dic = {}
    FR_wght_dic = {}
    
    FR_dic["Epitope Classes"] = Ab_One_Mut
    FR_wght_dic["Epitope Classes"] = Ab_One_Mut
    
    # compute mean IC50 per Ab_classes
    IC50_group = Escape_Fraction.groupby('condition', as_index=False).first()[['condition', 'IC50', 'group']]
    mean_IC50_per_group = IC50_group.groupby('group')['IC50'].mean().reset_index()
    print("Mean IC50 per Epitope Classes")
    print(mean_IC50_per_group)
    Mean_IC50 = np.ones(len(Ab_One_Mut))
    for i in range(len(Ab_One_Mut)):
        ab = Ab_One_Mut[i]
        if ab != "NTD":
            Mean_IC50[i] = (mean_IC50_per_group["IC50"].values[mean_IC50_per_group["group"] == ab])[0]
    
    for i in range(len(One_mut_lin_new)):
        if i != idx_WT:
            s = int(One_mut_lin_new[i])
            FR_dic[One_mut_lin_new[i]] = FR_Sites_Ab[:, i]
            if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                FR_wght_dic[One_mut_lin_new[i]] = FR_Sites_Ab[:, i]
            else:
                FR_wght_dic[One_mut_lin_new[i]] = FR_Sites_Ab[:, i]*Mean_IC50
    
    FR_df = pd.DataFrame(FR_dic)
    FR_df.to_csv(sys.argv[7])
    ### Save weighted FR DMS
    FR_df2 = pd.DataFrame(FR_wght_dic)
    FR_df2.to_csv(str(sys.argv[7])[:-4]+"_weighted.csv")
