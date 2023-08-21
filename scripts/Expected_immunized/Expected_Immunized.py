#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pickle
import re
import joblib as jb

from .util import Immunity_dynamics_fftconvolve, Antibody_ranges, Find_IC50_ranges, spikegroups_proportion

"""Load Infection Data"""
Population_Data_v0 = pd.read_csv('Data/caseAscertainmentTable.csv')
Population_Data = pd.read_csv("Data/caseAscertainmentTable_reportedCasesRatio.csv")

"""Load population data starting from July 1st, 2021"""
Population_Data = Population_Data.drop(index = Population_Data.index[:list(Population_Data['date']).index("2021-07-01")])

t = np.arange(1, len(Population_Data['date'])+1, 1) # array of timepoints at which to compute the antibody concentration
infection_data_corrected = Population_Data['minNTrue'].values
t_dates = Population_Data['date'].values
days_incidence = list(Population_Data['date']) 
Population_Data["pop"] = (Population_Data_v0["pop"].values[0])*np.ones(len(days_incidence))

file1 = open("Data/Cross_with_delta_validation.pck", "rb") # premade simulations
Cross_with_delta_validation = pickle.load(file1)
variant_x_names_show = Cross_with_delta_validation["variant_list"]
Cross_with_delta_validation.pop("variant_list")
file1.close()

file1 = open("Data/Cross_react_dic_spikegroups.pck", "rb") # Check that this is the file you want to load
Cross_react_dic = pickle.load(file1)
variant_x_names_cross = Cross_react_dic["variant_list"]
file1.close()

Ab_classes = list(Cross_react_dic.keys())

"""Compute Antibody concentration over time for a range of t_half and t_max"""
thalf_vec = np.linspace(25, 69, 15) 
tmax_vec = np.linspace(14, 28, 5)
c_t_vec, c_dframe_dic, dataname = Antibody_ranges(thalf_vec, tmax_vec, t, Ab_classes)
IC50xx_dic, mean_IC50xx_dic = Find_IC50_ranges(thalf_vec, tmax_vec, t, Ab_classes,  Cross_with_delta_validation)


file1 = open("Data/SpikeGroups.pck", "rb") 
SpikeGroups_list = pickle.load(file1)["names"]
file1.close()

def ei_util(i):
    variant_to_sim = [SpikeGroups_list[i]]
    EI = {}
    EI["Days"] = days_incidence
    
    for key in c_dframe_dic.keys():
        PK_dframe = c_dframe_dic[key]
        key_num = np.array(re.findall(r"\d+", key)).astype(int)

        Res_sub_0 = Immunity_dynamics_fftconvolve(t, PK_dframe, infection_data = infection_data_corrected, 
                                                     present_variant_list = SpikeGroups_list, 
                                                     tested_variant_list =  variant_to_sim, 
                                                     variant_name = variant_x_names_cross, 
                                                     variant_proportion =  spikegroups_proportion, 
                                                     Ab_classes = Ab_classes, 
                                                     IC50xx= mean_IC50xx_dic,
                                                     Cross_react_dic = Cross_react_dic, 
                                                     )
        
        EI["t_half = %.3f \nt_max = %.3f"%(thalf_vec[key_num[0]], tmax_vec[key_num[1]])] = Res_sub_0
        
    """ Save Dynamics Without Vaccination """
    EI_df = pd.DataFrame(EI)
    EI_df.to_csv("Data/Immunized_SpikeGroup_%s_all_PK.csv"%variant_to_sim[0])


try:
    jb.Parallel(n_jobs = -1)(jb.delayed(ei_util)(i) for i in range(len(SpikeGroups_list)))    
except:
    jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(ei_util)(i) for i in range(len(SpikeGroups_list)))    
