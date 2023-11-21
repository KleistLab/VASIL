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
    if joblib in (True, "joblib"):
        """Parallel codes --- macOS Monterey 12.5 crashes --- Not used by default """
        pfunc = partial(sub_Bind, tiled_esc = tiled_esc, Where_Mut = Where_Mut, Where_Cond = Where_Cond)
        status = False
        try:
            jb_res = list(jb.Parallel(n_jobs = -1, backend = "loky")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
            status = True
            print("run joblib.Parallel")
        except:
            try:
                jb_res = list(jb.Parallel(n_jobs = -1, backend = "multiprocessing")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
                status=True
                print("run joblib.Parallel")
            except:
                jb_res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
                status=True
                print("run joblib.Parallel")
        
        if status:
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
            mut_sites += list(np.array(mut_sites_per_variant[var]))
    
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



"""Compute Cross reactivity to Delta for valitation"""
# use Delta: B.1.617.2 for validation (Used in Clinical data paper)
variant_x_names_show = ["Wuhan-Hu-1", "Delta: B.1.617.2"]
mut_dic_show = {"Wuhan-Hu-1":[], "Delta: B.1.617.2": [614, 950, 142, 452, 681, 19, 478]}

Cross_with_delta_validation, Missed, Greater_one = cross_reactivity((variant_x_names_show,variant_x_names_show), 
															 Escape_Fraction, Ab_classes, 
                                                        	 mut_dic_show)

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
            for s in sites:
                s = int(s)
                if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                    FR_sites *= 10
            FR_NTD[i, j] = FR_sites
            FR_NTD[j, i] = FR_sites

Cross_with_delta_validation["NTD"] = FR_NTD
Cross_with_delta_validation["variant_list"] = variant_x_names_show
try:
    file0 = open(sys.argv[len(sys.argv)-2], "wb") 
    pickle.dump(Cross_with_delta_validation, file0)
    file0.close()
except:
    pass

"""Compute Cross reactivity"""
Cross_react_dic = {}
AB = ""
a = 1
print("Cross reactivity computation might take a while")
try:
    joblib = str(sys.argv[8])
except:
    joblib = None

"""Load lineage name to assess and it's mutation profile"""
try:
    n_groups = int(sys.argv[4]) 
    Lin_name = "Groups"
except:
    Lin_name = sys.argv[4]

if Lin_name not in ("ALL", "FR_DMS_sites", "missing"):
    if Lin_name != "Groups":
        try:
            mut_file = open(sys.argv[5], "r")
            mut_lin0 = mut_file.readlines()
            mut_file.close()

            mut_Lin = []
            for mut in mut_lin0:
                if mut[:3] not in ("DEL", "del"):
                    if len(re.findall(r'\d+', mut))>0:
                        mut_Lin.append(re.findall(r'\d+', mut)[0])       
                        mut_Lin = list(np.unique(np.array(mut_Lin).astype(str)))

            """Update mutation profile dictionary"""
            mut_x_sites_dic_updated = mut_x_sites_dic.copy()
            if Lin_name not in variant_x_names_cross: ### Keep Lin_name as it is
                mut_x_sites_dic_updated[Lin_name] = mut_Lin
            else:
                mut_x_sites_dic_updated[Lin_name] = mut_Lin
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
                               mut_x_sites_dic_updated)
                    
                    #Only the information for the specific lineage studied is required for immunological landscape calculation
                    #the FRxy_ab matrix is kept only for compatibility with other codes
                    
                    FRxy_ab[w_lin, g[s]] = Cross_Lin[ab][0, :]
                    FRxy_ab[g[s], w_lin] = Cross_Lin[ab][0, :]
                Cross_react_dic[ab] = FRxy_ab
            a +=1
        
        """Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""   
        print("Cross reactivity Epitope NTD")
        n = len(Cross_react_dic["variant_list"])
        FR_NTD = np.ones((n, n))
        for i in range(n):
            var_1 = Cross_react_dic["variant_list"][i]
            for j in range(n):
                if i > j:
                    var_2 = Cross_react_dic["variant_list"][j]
        
                    sites_1 = set(np.array(mut_x_sites_dic_updated[var_1]).astype(int))
                    sites_2 = set(np.array(mut_x_sites_dic_updated[var_2]).astype(int))
        
                    sites = list(sites_1.symmetric_difference(sites_2))
                    FR_sites = 1
                    for s in sites:
                        s = int(s)
                        if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                            FR_sites *= 10
                    FR_NTD[i, j] = FR_sites
                    FR_NTD[j, i] = FR_sites
        Cross_react_dic["NTD"] = FR_NTD
        file0 = open(sys.argv[7], "wb") 
        pickle.dump(Cross_react_dic, file0)
        file0.close()

    else: ### must be groups
        file = open("Spikegroups_membership.pck", "rb")
        Pseudogroup_dic = pickle.load(file)
        file.close()

        k = 5
        Lin_list = []
        while k<(n_groups+5):
            Lin_list.append(sys.argv[k])
            k +=1
        
        mut_sim = []
        while k<(n_groups+5+len(Lin_list)):
            mut_sim.append(sys.argv[k])
            k +=1
        
        mut_x_sites_dic_updated = mut_x_sites_dic.copy()
        for j in range(len(Lin_list)):
            if Lin_list[j] not in list(Pseudogroup_dic.keys()):
                mut_file = open(mut_sim[j], "r")
                mut_lin0 = mut_file.readlines()
                mut_file.close()
                mut_Lin = []
                for mut in mut_lin0:
                    if mut[:3] not in ("DEL", "del"):
                        if len(re.findall(r'\d+', mut))>0:
                            mut_Lin.append(re.findall(r'\d+', mut)[0])       
                            mut_Lin = list(np.unique(np.array(mut_Lin).astype(str)))
                """Update mutation profile dictionary"""
                mut_x_sites_dic_updated[Lin_list[j]] = mut_Lin
            else:
                mut_x_sites_dic_updated[Lin_list[j]] = mut_x_sites_dic[Pseudogroup_dic[Lin_list[j]]]
                
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
        
        status_sim = []
        for i in range(len(Lin_list)):
            Cross_i = {}
            Cross_i["variant_list"] = list(variant_x_names_cross)+ [Lin_list[i]]
            
            if Lin_list[i] not in list(variant_x_names_cross):
                w_lin = len(variant_x_names_cross)
                Cross_i["variant_list"] = list(variant_x_names_cross) + [Lin_list[i]]
            else:
                w_lin = list(variant_x_names_cross).index(Lin_list[i])
                Cross_i["variant_list"] = list(variant_x_names_cross)
            
            try:
                file_test = open("results/Cross_react_dic_spikegroups_ALL.pck", "rb")
                file_test.close()
                extract = True 
            except:
                extract = False # file is not present and thus if must be recomputed
            
            if (Lin_list[i] not in list(Pseudogroup_dic.keys())) or (not extract):
                a = 1
                for ab in Ab_classes:  
                    if ab!= "NTD":
                        if Lin_list[i] not in list(variant_x_names_cross):
                            FRxy_ab = np.ones((len(variant_x_names_cross)+1, len(variant_x_names_cross)+1))
                        else:
                            FRxy_ab = np.ones((len(variant_x_names_cross), len(variant_x_names_cross)))
                            
                        FRxy_ab = np.ones((len(variant_x_names_cross)+1, len(variant_x_names_cross)+1))
                        print("Assess lineage %s| %d out of %d with the NTD-RBD mutation positions"%(Lin_list[i], i+1,len(Lin_list)), mut_x_sites_dic_updated[Lin_list[i]])
                        print("Cross reactivity Epitope %s, countdown"%ab, a, "out of %d epitope clases"%len(Ab_classes)) 
                        for s in range(len(g)):
                            Cross_Lin, Missed, Greater_one = cross_reactivity(([Lin_list[i]], g_var[s]), 
                                       Escape_Fraction, 
                                       [ab],
                                       mut_x_sites_dic_updated)
                            
                           
                            #Only the information for the specific lineage studied is required for immunological landscape calculation
                            #the FRxy_ab matrix is kept only for compatibility with other codes
                            
                            FRxy_ab[w_lin, g[s]] = Cross_Lin[ab][0, :]
                            FRxy_ab[g[s], w_lin] = Cross_Lin[ab][0, :]
                
                        Cross_i[ab] = FRxy_ab
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
                
                            sites_1 = set(np.array(mut_x_sites_dic_updated[var_1]).astype(int))
                            sites_2 = set(np.array(mut_x_sites_dic_updated[var_2]).astype(int))
                
                            sites = list(sites_1.symmetric_difference(sites_2))
                            FR_sites = 1
                            for s in sites:
                                s = int(s)
                                if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                                    FR_sites *= 10
                            FR_NTD[i1, j1] = FR_sites
                            FR_NTD[j1, i1] = FR_sites
                Cross_i["NTD"] = FR_NTD
            else:
                ### open global cross_reactivity file which must be present
                a = 1  
                print("Assess lineage %s| %d out of %d with the NTD-RBD mutation positions"%(Lin_list[i], i+1,len(Lin_list)), mut_x_sites_dic_updated[Lin_list[i]])
                print("Load : %s is present in general file results/Cross_react_dic_spikegroups_ALL.pck"%Lin_list[i]) 
                file_c = open("results/Cross_react_dic_spikegroups_ALL.pck", "rb") 
                Cross_global = pickle.load(file_c)
                variant_global = Cross_global["variant_list"]
                Cross_global.pop("variant_list")
                Ab_global = Cross_global.keys()
                file_c.close()
                if Lin_list[i] not in list(variant_x_names_cross):
                    w_lin = len(variant_x_names_cross)
                    Cross_i["variant_list"] = list(variant_x_names_cross) + [Lin_list[i]]
                else:
                    w_lin = list(variant_x_names_cross).index(Lin_list[i])
                    Cross_i["variant_list"] = list(variant_x_names_cross)
                    
                for ab in Ab_global:  
                    if Lin_list[i] not in list(variant_x_names_cross):
                        FRxy_ab = np.ones((len(variant_x_names_cross)+1, len(variant_x_names_cross)+1))
                    else:
                        FRxy_ab = np.ones((len(variant_x_names_cross), len(variant_x_names_cross)))
                        
                    for u1 in range(len(variant_x_names_cross)):
                        v_u1 = variant_x_names_cross[u1]
                        if v_u1 in variant_global:
                            FRxy_ab[w_lin, u1] = Cross_global[ab][list(variant_global).index(Lin_list[i]), list(variant_global).index(v_u1)]
                            FRxy_ab[u1, w_lin] = FRxy_ab[w_lin, u1]
                        else:
                            # recompute it, should not happen normaly
                            # usefull when using previously computed cross reactivity file where all covsonar lineages (not only spikegroups) are present (assigned the FR of their spikegroups)
                            if ab != "NTD":
                                Cross_Lin, Missed, Greater_one = cross_reactivity(([Lin_list[i]], [v_u1]), 
                                           Escape_Fraction, 
                                           [ab],
                                           mut_x_sites_dic_updated)
                                
                               
                                #Only the information for the specific lineage studied is required for immunological landscape calculation
                                #the FRxy_ab matrix is kept only for compatibility with other codes
                                FRxy_ab[w_lin, u1] = Cross_Lin[ab][0, 0]
                                FRxy_ab[u1, w_lin] = FRxy_ab[w_lin, u1]
                            else:
                                sites_1 = set(np.array(mut_x_sites_dic_updated[Lin_list[i]]).astype(int))
                                sites_2 = set(np.array(mut_x_sites_dic_updated[v_u1]).astype(int))
                                sites = list(sites_1.symmetric_difference(sites_2))
                                FR_sites = 1
                                for s in sites:
                                    s = int(s)
                                    if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                                        FR_sites *= 10
                                FRxy_ab[w_lin, u1] = FR_sites
                                FRxy_ab[u1, w_lin] = FR_sites
                    Cross_i[ab] = FRxy_ab
                    a +=1 
                
            status_sim.append("Done")
            file0 = open(sys.argv[k+1]+"/Cross_%s.pck"%Lin_list[i], "wb") 
            pickle.dump(Cross_i, file0)
            file0.close()
        stat_df = pd.DataFrame({"Lineages":Lin_list, "computed_cross":status_sim})
        stat_df.to_csv(sys.argv[k+1]+"/cross_status.csv")

elif Lin_name == "ALL":  
    """Break runs into manageable pieces"""
    g = []
    g_var =[]
    inds = np.arange(0, len(variant_x_names_cross)).astype(int)
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
    
          
    for ab in Ab_classes:
        if ab!= "NTD":
            FRxy_ab = np.ones((len(variant_x_names_cross), len(variant_x_names_cross)))
            num_to_do = 0
            for s1 in range(len(g)):
                print("Assess all spikegroups with the NTD-RBD mutation positions ")
                print("Cross reactivity Epitope %s, countdown"%ab, a, "out of %d epitope clases"%len(Ab_classes))
                print("Run Cross for %d out of %d (to achieve %d/%d spikegroups)"%(s1+1, len(g), num_to_do+len(g[s1]), len(variant_x_names_cross)))
                sub_FR = np.ones((len(g[s1]), len(variant_x_names_cross)))
                for s2 in range(len(g)):
                    Cross_Lin, Missed, Greater_one = cross_reactivity((g_var[s1], g_var[s2]), 
                                                                       Escape_Fraction, 
                                                                       [ab],
                                                                       mut_x_sites_dic, joblib=True)

                    sub_FR[:, g[s2]] = Cross_Lin[ab]
        
                FRxy_ab[g[s1], :] = sub_FR
                num_to_do +=len(g[s2])
            Cross_react_dic[ab] = FRxy_ab
        a +=1
    
    
    Cross_react_dic["variant_list"] = list(variant_x_names_cross)
    """Add FR to NTD-targeting AB assuming a FR of 10 to each mutations sites included in NTD Antigenic supersite"""  
    print("Cross reactivity spikegroups for Epitope NTD")
    n = len(Cross_react_dic["variant_list"])
    FR_NTD = np.ones((n, n))
    for i in range(n):
        var_1 = Cross_react_dic["variant_list"][i]
        for j in range(n):
            if i > j:
                var_2 = Cross_react_dic["variant_list"][j]
    
                sites_1 = set(np.array(mut_x_sites_dic[var_1]).astype(int))
                sites_2 = set(np.array(mut_x_sites_dic[var_2]).astype(int))
    
                sites = list(sites_1.symmetric_difference(sites_2))
                FR_sites = 1
                for s in sites:
                    s = int(s)
                    if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                        FR_sites *= 10
                FR_NTD[i, j] = FR_sites
                FR_NTD[j, i] = FR_sites
        
    Cross_react_dic["NTD"] = FR_NTD
    file0 = open(sys.argv[7], "wb") 
    pickle.dump(Cross_react_dic, file0)
    file0.close()

elif Lin_name == "missing":
    file_c = open("results/Cross_react_dic_spikegroups_ALL.pck", "rb") 
    Cross_global = pickle.load(file_c)
    variant_global = list(Cross_global["variant_list"])
    Cross_global.pop("variant_list")
    Ab_global = Cross_global.keys()
    file_c.close()
    
    """Find indexes of missing and not missing variants"""
    Lin_miss = []
    loc_not_miss = []
    loc_in_cross = []
    for lin in variant_x_names_cross:
        if lin not in variant_global:
            Lin_miss.append(lin)
        else:
            loc_in_cross.append(list(variant_x_names_cross).index(lin))
            loc_not_miss.append(list(variant_global).index(lin))
            
    if len(Lin_miss) == 0:
        Cross_react_dic = Cross_global.copy()
    else:
        w_in_cross = np.arange(0, len(variant_x_names_cross)).astype(int)[np.array(loc_in_cross)]
        variants_in_global = np.array(variant_x_names_cross)[w_in_cross]
        Cross_react_dic["variant_list"] = list(variants_in_global) + Lin_miss
        a = 1
        g = []
        g_var =[]
        inds = np.arange(0, len(variants_in_global)).astype(int)
        cut_step = 100
        if len(variants_in_global)>cut_step:
            cut1 = 0
            cut2 = cut_step
            while cut2<len(variants_in_global):
                g.append(inds[cut1:cut2])
                g_var.append(list(np.array(variants_in_global)[cut1:cut2]))
                cut1=cut2
                cut2+=min(cut_step, len(variants_in_global)-cut2)
            g.append(inds[cut1:cut2])
            g_var.append(list(np.array(variants_in_global)[cut1:cut2]))
        else:
            g.append(inds)
            g_var.append(variants_in_global)

        w_global = np.arange(0, len(variant_global)).astype(int)[np.array(loc_not_miss)] ## location of variant_x_names cross in global file
        for ab in Ab_global:
            Cross_react_dic[ab] = np.ones((len(variants_in_global)+len(Lin_miss), len(variants_in_global)+len(Lin_miss)))
            if ab != "NTD":
                w_miss = len(variants_in_global)
                for s in range(len(g)):
                    print("Assess missing | num %d vs (%d, %d (max %d)) with the NTD-RBD mutation positions"%(len(Lin_miss), s*cut_step, min((s+1)*cut_step, len(variants_in_global)),len(variants_in_global)))
                    print("Cross reactivity Epitope %s, countdown"%ab, a, "out of %d epitope clases"%len(Ab_global)) 
                    Cross_Lin, Missed, Greater_one = cross_reactivity((Lin_miss, g_var[s]), 
                               Escape_Fraction, 
                               [ab],
                               mut_x_sites_dic,
                               joblib = True)
                    
                    Cross_react_dic[ab][w_miss:, :w_miss][:, g[s]] = Cross_Lin[ab]
                    Cross_react_dic[ab][:w_miss, w_miss:][g[s], :] = Cross_Lin[ab].T
                
                print("Assess %d missing vs. %d missing with the NTD-RBD mutation positions"%(len(Lin_miss), len(Lin_miss)))
                Cross_Lin, Missed, Greater_one = cross_reactivity((Lin_miss, Lin_miss), 
                           Escape_Fraction, 
                           [ab],
                           mut_x_sites_dic,
                           joblib = True)
                
                Cross_react_dic[ab][w_miss:, w_miss:] = Cross_Lin[ab]
            else:
                n = len(Cross_react_dic["variant_list"])
                FR_NTD = np.ones((n, n))
                for i in range(len(variants_in_global), n):
                    var_1 = Cross_react_dic["variant_list"][i]
                    for j in range(n):
                        var_2 = Cross_react_dic["variant_list"][j]
            
                        sites_1 = set(np.array(mut_x_sites_dic[var_1]).astype(int))
                        sites_2 = set(np.array(mut_x_sites_dic[var_2]).astype(int))
            
                        sites = list(sites_1.symmetric_difference(sites_2))
                        FR_sites = 1
                        for s in sites:
                            s = int(s)
                            if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)):
                                FR_sites *= 10
                        FR_NTD[i, j] = FR_sites
                        FR_NTD[j, i] = FR_sites
                    
                Cross_react_dic[ab] = FR_NTD
            
            Cross_react_dic[ab][:len(variants_in_global), :len(variants_in_global)] = Cross_global[ab][w_global, :][:, w_global]
            a +=1 
    
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
                                            joblib=joblib)
        
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
            for s in sites:
                s = int(s)
                if ((14<=s)&(s<=20)) or ((140<=s)&(s<=158)) or ((245<=s)&(s<=264)): ### Antigenic supersites
                    FR_sites *= 10
            FR_NTD[i] = FR_sites
    
    FR_Sites_Ab = np.row_stack((FR_Sites_Ab, FR_NTD)) 
    
    
    ### Saving file
    FR_dic = {}
    FR_dic["Epitope Classes"] = Ab_One_Mut
    for i in range(len(One_mut_lin_new)):
        if i != idx_WT:
            FR_dic[One_mut_lin_new[i]] = FR_Sites_Ab[:, i]
            
    FR_df = pd.DataFrame(FR_dic)
    FR_df.to_csv(sys.argv[7]) 

    
    