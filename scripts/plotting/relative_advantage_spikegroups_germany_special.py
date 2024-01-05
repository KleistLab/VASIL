import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import pdb
import seaborn as sns
import pickle
import mpl_axes_aligner
import re
import scipy.stats

def moving_average(X, window = 7):
    u = np.zeros(len(X))
    u[:window] = X[:window]
    for i in range(window, len(X)):
        u[i] = np.mean(X[i-window:i+1])
    
    return u

def lower_ci(z_crit, proportion, total_size):
    B2 = total_size
    B6 = proportion
    B7 = z_crit
    B10 = z_crit**2/total_size
    return((B6+B10/2)/(B10+1)-B7/(B10+1)*np.sqrt(B6*(1-B6)/B2+B10/(4*B2)))

def upper_ci(z_crit, proportion, total_size):
    B2 = total_size
    B6 = proportion
    B7 = z_crit
    B10 = z_crit**2/total_size
    return((B6+B10/2)/(B10+1)+B7/(B10+1)*np.sqrt(B6*(1-B6)/B2+B10/(4*B2)))

ES_lin_dir = sys.argv[1]
S_mean_file = sys.argv[2]

lineage_freq = pd.read_csv(sys.argv[3])
threshold = float(sys.argv[4])

# needs to be updated to allow individual weighting 
S_mean_df = pd.read_csv(S_mean_file)
S_all_mean = S_mean_df.to_numpy()[:, (S_mean_df.columns != "Days")&(S_mean_df.columns != "Unnamed: 0")].astype(float)
t_dates = S_mean_df["Days"].tolist()

# processing of frequency data
try:
    lineage_freq.drop(columns = "Unnamed: 0", inplace = True)
except:
    pass

### Drop proportions column that does not start with correspond to E[Susceptible] files
lineage_freq = lineage_freq.drop(index = lineage_freq.index[:list(lineage_freq['date']).index(S_mean_df["Days"][0])])
freqs = lineage_freq.loc[:, lineage_freq.columns != 'date']
# imputing frequencies below threshold and normalization
prop_mask = np.all(lineage_freq.loc[:, lineage_freq.columns != 'date'] == 0.0, axis = 1)
freqs = freqs.mask(freqs < threshold)
freqs = freqs.fillna(0)
col_sums = freqs.sum(axis = 1).values
freqs = freqs.divide(col_sums, axis="rows")
freqs = freqs.fillna(0)
lineage_freq.loc[:, lineage_freq.columns != 'date'] = freqs
day_prop = lineage_freq["date"].tolist()
t_prop = np.arange(len(day_prop)).astype(int)

# proportion for the plot after end of observation
t_after = np.array([0, 30, 60, 90]) + len(t_dates) 
lineage_freq_after = lineage_freq[lineage_freq["date"]> "2023-04-15"]
lineage_freq_after["month"] = ((lineage_freq_after["date"].str[5:7]).astype("int"))
lineage_freq_after = lineage_freq_after.drop("date", axis = 1)

lineage_freq_after_monthly = lineage_freq_after.groupby(by=["month"]).sum()
month_sum = np.sum(lineage_freq_after_monthly, axis = 1)
lineage_freq_after_monthly = lineage_freq_after_monthly.reset_index()

# for confidence intervals
alpha = 0.05
z_crit = scipy.stats.norm.ppf(1-alpha/2)

# import Stichprobe for genome count
stichprobe = pd.read_csv("Stichprobe_RKI-Jul2021_to_9thAug2023_KW.tsv", sep='\t')
stichprobe["month"] = ((stichprobe["date"].str[5:7]).astype("int"))
stichprobe = stichprobe[stichprobe["date"] > "2023-04-01"]
genomes_month = stichprobe.groupby(by=["month"])["month"].count()

import matplotlib
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)

file = open("Spikegroups_membership.pck", "rb")
Pseudogroup_dic = pickle.load(file)
file.close()    
def plot_fit(ES_df_dir, lineage_list, color_list, w_save = len(sys.argv)-1, already_prop = np.zeros((len(t_prop)))):
    # plotting
    PreFig(xsize = 20, ysize = 20)
    fig = plt.figure(figsize = (15, 7))
    ax = fig.add_subplot(1, 1, 1)
    ### end of observation line
    ax.axvline(x = len(t_dates) - 1, ymin = -1, ymax = 1, ls = "--", linewidth = 2, color = "grey")
    # different axis for proportions
    ax_twin = ax.twinx()
    
    ### Separate figure for proportions
    PreFig(xsize = 20, ysize = 20)
    fig_prop = plt.figure(figsize = (9, 7))
    ax_prop = fig_prop.add_subplot(1, 1, 1)
    
    status_list = []
    Pseudo_done_global = []
    
    for k in range(len(lineage_list)):
        if lineage_list[k][-4:] == ".ALL" and not re.search("/", lineage_list[k]):
            splited_var = []
            for x in list(Pseudogroup_dic.keys()):
                if (x[:len(lineage_list[k][:-4])] == lineage_list[k][:-4]):
                    if len(x) == len(lineage_list[k][:-4]):
                        splited_var.append(x)
                    elif len(x)>len(lineage_list[k][:-4]):
                        if x[len(lineage_list[k][:-4])] == ".":
                            splited_var.append(x)
                    
        else:
            splited_var = np.array(lineage_list[k].split("/"))
            splited_var = splited_var[~(splited_var == "")]
            splited_var = splited_var[~(splited_var == " ")]
         
        splited_var0 = np.array(splited_var).copy()
        splited_var = []
        for var in splited_var0:
            if var[-4:] == ".ALL":
                for x in list(Pseudogroup_dic.keys()):
                    if x[:len(var[:-4])] == var[:-4]:
                        if len(x) == len(var[:-4]):
                            splited_var.append(x)
                        elif len(x)>len(var[:-4]):
                            if x[len(var[:-4])] == ".":
                                splited_var.append(x)
            else:
                splited_var.append(var)
        
        num_avail = 0
        # plotting
        PreFig(xsize = 20, ysize = 20)
        fig_k = plt.figure(figsize = (9, 7))
        ax_k = fig_k.add_subplot(1, 1, 1)
        ### end of observation line
        ax_k.axvline(x = len(t_dates) - 1, ymin = -1, ymax = 1, ls = "--", linewidth = 2, color = "grey")

        # different axis for proportions
        ax_k_twin = ax_k.twinx()
        
        
        lab_done = {}
        lab_done[lineage_list[k]] = ""
        Pseudo_done = {}
        Pseudo_done[lineage_list[k]] = ""
        plot_prop = []
        Pseudo_Prop = np.zeros((len(t_prop)))
        Prop_after = np.zeros((len(month_sum.values)))
        start = 0
        ES_list = []
        prop_list = []
        for x in range(len(splited_var)):
            lineage = splited_var[x]                
            try:
                if Pseudogroup_dic[lineage] not in (Pseudo_done[lineage_list[k]].split("/")):
                    ES_df = pd.read_csv(ES_df_dir+"/Susceptible_SpikeGroup_%s_all_PK.csv"%Pseudogroup_dic[lineage])
                    num_avail +=1 
                    
                    try:
                        ES_df.drop(columns = "Unnamed: 0", inplace = True)
                    except:
                        pass
                    
                    es_cols = ES_df.columns
                    ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                    ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))    
                    run = True
            except:
                try:
                    if "Placeholder"+ lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                        ES_df = pd.read_csv(ES_df_dir+"/Susceptible_SpikeGroup_%s_all_PK.csv"%lineage)
                        num_avail +=1
                        try:
                            ES_df.drop(columns = "Unnamed: 0", inplace = True)
                        except:
                            pass
                        
                        es_cols = ES_df.columns
                        ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                        ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))
                        run = True
                except:
                    try:
                        if Pseudogroup_dic[lineage] not in (Pseudo_done[lineage_list[k]].split("/")):
                            ES_df = pd.read_csv("results/Immunological_Landscape_ALL/Susceptible_SpikeGroup_%s_all_PK.csv"%Pseudogroup_dic[lineage])
                            num_avail +=1
                            try:
                                ES_df.drop(columns = "Unnamed: 0", inplace = True)
                            except:
                                pass
                            
                            es_cols = ES_df.columns
                            ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                            ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))
                            run = True
                    except:
                        try:
                            if Pseudogroup_dic[lineage] not in (Pseudo_done[lineage_list[k]].split("/")):
                                ES_df = pd.read_csv("results/Immunological_Landscape/Susceptible_SpikeGroup_%s_all_PK.csv"%Pseudogroup_dic[lineage])
                                num_avail +=1
                            
                                try:
                                    ES_df.drop(columns = "Unnamed: 0", inplace = True)
                                except:
                                    pass
                                
                                es_cols = ES_df.columns
                                ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                                ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))
                                run = True
                        except:
                             try:
                                 if "Placeholder"+ lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                                     ES_df = pd.read_csv("results/Immunological_Landscape/Susceptible_SpikeGroup_%s_all_PK.csv"%lineage)
                                     num_avail +=1
                                     try:
                                         ES_df.drop(columns = "Unnamed: 0", inplace = True)
                                     except:
                                         pass
                                        
                                     es_cols = ES_df.columns
                                     ES_df = ES_df[ES_df['Days'].isin(t_dates)]
                                     ES_list.append(ES_df.to_numpy()[:, es_cols!="Days"].astype(float))
                                    
                                     run = True
                             except:
                                print("Computation needed: Expected Susceptible file is not available for %s"%lineage)
                                run = False
            
            # processing of Proportions data
            if run:   
                # change in relative frequency from genomic surveillance data 
                if lineage in list(Pseudogroup_dic.keys()):
                    lab_k = lineage + "*"+"/"
                    plot_prop.append(True)
                    if "Spike. " + Pseudogroup_dic[lineage] in lineage_freq.columns.astype(str):
                        if Pseudogroup_dic[lineage] not in (Pseudo_done[lineage_list[k]].split("/")):
                            sub_prop = moving_average(lineage_freq["Spike. " + Pseudogroup_dic[lineage]], window = 14)
                            sub_prop_after = lineage_freq_after_monthly["Spike. " + Pseudogroup_dic[lineage]]/month_sum.values
                            prop_list.append(sub_prop)
                            Pseudo_Prop += sub_prop
                            Prop_after += sub_prop_after
                            lab_done[lineage_list[k]] = lab_done[lineage_list[k]][:-1] + " + "+lab_k
                            if Pseudogroup_dic[lineage] not in Pseudo_done_global: ### Repeating Pseudogroups in variants combinations cannot be accounted twice in final proportion plot
                                already_prop += sub_prop
                                Pseudo_done_global.append(Pseudogroup_dic[lineage])
                        else:
                            lab_done[lineage_list[k]] += lab_k
                    elif Pseudogroup_dic[lineage] in lineage_freq.columns.astype(str):
                        if Pseudogroup_dic[lineage] not in (Pseudo_done[lineage_list[k]].split("/")):
                            sub_prop = moving_average(lineage_freq[Pseudogroup_dic[lineage]], window = 14)
                            sub_prop_after = lineage_freq_after_monthly[Pseudogroup_dic[lineage]]/month_sum.values
                            prop_list.append(sub_prop)
                            Pseudo_Prop += sub_prop
                            Prop_after += sub_prop_after
                            lab_done[lineage_list[k]] = lab_done[lineage_list[k]][:-1] + " + "+lab_k
                            if Pseudogroup_dic[lineage] not in Pseudo_done_global: ### Repeating Pseudogroups in variants combinations cannot be accounted twice in final proportion plot
                                already_prop += sub_prop
                                Pseudo_done_global.append(Pseudogroup_dic[lineage])
                        else:
                            lab_done[lineage_list[k]] += lab_k
                    
                    if Pseudogroup_dic[lineage] not in (Pseudo_done[lineage_list[k]].split("/")):
                        Pseudo_done[lineage_list[k]] += Pseudogroup_dic[lineage] + "/"
                    
                else:
                    lab_k = lineage+"/"
                    if lineage in lineage_freq.columns.astype(str):
                        plot_prop.append(True)
                        if "Placeholder"+lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                            sub_prop = moving_average(lineage_freq[lineage], window = 14)
                            sub_prop_after = lineage_freq_after_monthly[lineage]/month_sum.values
                            prop_list.append(sub_prop)
                            Pseudo_Prop += sub_prop
                            Prop_after += sub_prop_after
                            lab_done[lineage_list[k]] = lab_done[lineage_list[k]][:-1] + " + "+lab_k
                            if Pseudogroup_dic[lineage] not in Pseudo_done_global: ### Repeating Pseudogroups in variants combinations cannot be accounted twice in final proportion plot
                                already_prop += sub_prop
                                Pseudo_done_global.append(lineage)
                        else:
                            lab_done[lineage_list[k]] += lab_k

                    else:
                        plot_prop.append(False)
                        Pseudo_Prop += ma.masked_array(np.zeros(len(t_prop)), mask = np.ones(len(t_prop), dtype = bool))
                        if "Placeholder"+lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                            lab_done[lineage_list[k]] += " + "+lab_k
                            prop_list.append(ma.masked_array(np.zeros(len(t_prop)), mask = np.ones(len(t_prop), dtype = bool)))
                        else:
                            lab_done[lineage_list[k]] += " + "+lab_k
                    
                    if "Placeholder"+lineage not in (Pseudo_done[lineage_list[k]].split("/")):
                        Pseudo_done[lineage_list[k]]+="Placeholder"+ lineage +"/"
                
                start +=1
        
        lab_status = lab_done[lineage_list[k]][:-1]
        if not re.search(".ALL",lineage_list[k]):
            lab_k = lab_done[lineage_list[k]][:-1]
        else:
            lab_k = lineage_list[k]
        
        if lab_k[:3] == " + ":
            lab_k = lab_k[3:]
            
        if lab_status[:3] == " + ":
            lab_status = lab_status[3:]    
            
        lab_k_fn = (lab_k.replace("/", "_")).replace("*","").replace("+", "_") ## for filename 
        if len(lab_k_fn) > 10: # can't be too long
            lab_k_fn = lab_k_fn[:10] + "_et_al"
         
            
        if (lab_k != "" and num_avail != 0):
            if (num_avail == len(ES_list)):
                #ES_ranges= np.mean(np.array(ES_list), axis = 0) # compute the mean
                Props = np.array(prop_list)[:, :ES_list[0].shape[0]]
                Props_normed = np.divide(Props, np.sum(Props, axis = 0)[np.newaxis, :], out = np.zeros(Props.shape), where = np.sum(Props, axis = 0)[np.newaxis, :]!=0)
                ES_ranges = np.sum(np.array(ES_list)*Props_normed[:, :, np.newaxis], axis = 0) ## weighted mean
                ES_ranges[np.sum(Props, axis = 0) == 0, :] = np.mean(np.array(ES_list), axis = 0)[np.sum(Props, axis = 0) == 0, :]
            else:
                sys.exit("Loaded E[Susceptible] files were more than what is available for groups %s, recheck the loading process script/plotting/relative_advantage_groups.py Line 122-196"%lineage_list[k])
            
            # calculation of change in relative frequency from model
            Pseudo_Prop2 = Pseudo_Prop.copy()
            Pseudo_Prop2[Pseudo_Prop2 < 0.05] = 0
            Pseudo_Prop2 = list(Pseudo_Prop2)
            Pseudo_Prop_aligned = np.zeros(len(t_dates))
            already_prop_aligned = np.zeros(len(t_dates))
            day_prop_aligned = []
            t_prop_aligned = np.zeros(len(t_dates))
            prop_mask_aligned = np.zeros(len(t_dates)).astype(bool)
            gamma_prop = np.zeros(len(t_dates))
            SI_mask = np.zeros(len(t_dates)).astype(bool)
            for l in range(len(t_dates)):
                #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
                ## Dates are already set to always be to be successive
                if t_dates[l] in day_prop:
                    w_l = list(day_prop).index(t_dates[l])
                    Pseudo_Prop_aligned[l] = Pseudo_Prop[w_l]
                    prop_mask_aligned[l] = prop_mask.tolist()[w_l]
                    already_prop_aligned[l] = already_prop[w_l]
                    day_prop_aligned.append(t_dates[l])
                    t_prop_aligned[l] = w_l
                    try:
                        if Pseudo_Prop2[w_l] == 0 or Pseudo_Prop2[w_l+1] == 0:
                            gamma_prop[l] = float('nan')
                            SI_mask[l] == True
                        else:
                            gamma_prop[l] = (Pseudo_Prop2[w_l+1]/Pseudo_Prop2[w_l]) - 1
                    except:
                        gamma_prop[l] = float('nan')
                        SI_mask[l] = True
                else:
                    gamma_prop[l] = float('nan')
                    SI_mask[l] = True
            
            ### Making sure to aligne proportions and incidence timelines
            if len(day_prop_aligned) < len(day_prop):
                miss_days = [day for day in day_prop if day not in day_prop_aligned]
                day_prop_aligned = list(day_prop_aligned) + [day_prop[list(day_prop).index(miss_days[i])] for i in range(len(miss_days))]
                prop_mask_aligned = list(prop_mask_aligned) + [prop_mask.tolist()[list(day_prop).index(miss_days[i])] for i in range(len(miss_days))]
                t_prop_aligned = list(t_prop_aligned) + list(len(t_prop_aligned) + np.arange(len(miss_days)))
                Pseudo_Prop_aligned = list(Pseudo_Prop_aligned) + [Pseudo_Prop[list(day_prop).index(miss_days[i])] for i in range(len(miss_days))]
                already_prop_aligned = list(already_prop_aligned) + [already_prop[list(day_prop).index(miss_days[i])] for i in range(len(miss_days))]
                
                day_prop_aligned = list(day_prop_aligned)
                prop_mask_aligned = np.array(prop_mask_aligned)
                t_prop_aligned = np.array(t_prop_aligned)
                Pseudo_Prop_aligned = np.array(Pseudo_Prop_aligned)
                already_prop_aligned = np.array(already_prop_aligned)
             
            # calculation of relative fitness
            gamma_SI = np.zeros((len(t_dates), ES_ranges.shape[1]))
            for i in range(ES_ranges.shape[1]):
                S_x = ES_ranges[:, i]
                S_mean = S_all_mean[:, i]
                gamma_SI[:, i] = np.divide(S_x - S_mean, S_mean, out = S_x, where = (S_mean != 0)&(S_x!=0))
            
            
            # get min max gamma over PK at each timepoints
            inds_dates = np.arange(0,len(t_dates),1)
            SI_mask = np.array(SI_mask) + prop_mask_aligned[:len(inds_dates)] ### props are already aligned with indicence date
            gamma_SI_min, gamma_SI_max = np.min(gamma_SI, axis = 1), np.max(gamma_SI, axis = 1)
            gamma_SI_max = ma.masked_array(gamma_SI_max, mask = SI_mask)
            gamma_SI_min = ma.masked_array(gamma_SI_min, mask = SI_mask)
            Pseudo_Prop_masked = ma.masked_array(Pseudo_Prop, mask = prop_mask)

            # calculate CIs
            try:
                Prop_after = Prop_after.to_numpy()
            except:
                1+1
                
            Prop_after_lower = Prop_after - lower_ci(z_crit, Prop_after, genomes_month.values)
            Prop_after_lower[Prop_after_lower <0 ] = 0
            Prop_after_upper = Prop_after + upper_ci(z_crit, Prop_after, genomes_month.values)
            asymmetric_error = [Prop_after_lower, Prop_after_lower]
    
            x_diff = (-1)**k * 0.5
    
            ax_twin.errorbar(t_after + x_diff, Prop_after, yerr=asymmetric_error, fmt='o', color = color_list[k])
            ax_twin.plot(t_after, Prop_after, color = color_list[k], label = lab_k)
    

            ax.fill_between(inds_dates, gamma_SI_min, gamma_SI_max, color = color_list[k], alpha = 0.2)
            ax.plot(inds_dates, (gamma_SI_min + gamma_SI_max)/2, color = color_list[k], label = lab_k, linewidth = 3)
            #ax_twin.plot(t_prop_aligned, 100*Pseudo_Prop_masked_group, linewidth = 3, color = color_list[k], label = lab_k)
            #ax_twin.plot(t_after, 100*Prop_after, linewidth = 3, color = color_list[k], label = lab_k)
            #ax_twin.scatter(t_prop_aligned, Pseudo_Prop_masked, marker = ".", color = color_list[k])

            ### Plot spikegroups frequencies
            if np.all(plot_prop):
                Pseudo_Prop_masked = ma.masked_array(Pseudo_Prop_aligned, mask=prop_mask_aligned)
                ax_prop.plot(t_prop_aligned, 100*Pseudo_Prop_masked, linewidth = 3, color = color_list[k], label = lab_k)
                #ax_prop.scatter(t_prop_aligned, 100*Pseudo_Prop_masked, marker = ".", color = color_list[k])

            
            ax_k.fill_between(inds_dates, gamma_SI_min, gamma_SI_max, color = color_list[k], alpha = 0.3, label = lab_k) 
            ax_k_twin.plot(t_prop_aligned, 100*Pseudo_Prop_masked, linewidth = 3, color = color_list[k], label = lab_k) 
            #ax_k_twin.scatter(t_prop_aligned, Pseudo_Prop_masked, marker = ".", color = color_list[k])
            ax_k.axhline(xmin = 0, xmax = len(t_dates), ls = "--", linewidth = 2, color = "black")
            
            ymin1, ymax1 = ax_k.get_ylim()
            ymin2, ymax2 = ax_k_twin.get_ylim()
            #ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            ax_k.set_ylim((ymin1, ymax1)) 
            #ax_k_twin.set_ylim((ymin, ymax))
            #loc0k = min(np.abs(ymin1)/(np.abs(ymin1)+np.abs(ymax1)), np.abs(ymax1)/(np.abs(ymin1)+np.abs(ymax1)))
            #mpl_axes_aligner.align.yaxes(ax_k, 0, ax_k_twin, 0, loc0k)
            
            try:
                if str(sys.argv[5]) in list(t_dates):
                    x_min = list(t_dates).index(str(sys.argv[5]))
                else:
                    x_min = 0
                if str(sys.argv[5]) in list(day_prop_aligned):
                    x_min1 = day_prop_aligned.index(str(sys.argv[5]))
                else:
                    x_min1 = 0
                if str(sys.argv[6]) in list(day_prop_aligned):
                    x_max1 = day_prop_aligned.index(str(sys.argv[6]))
                else:
                    x_max1 = len(day_prop_aligned) - 1
                if str(sys.argv[6]) not in list(t_dates):
                    x_max = (len(t_dates) - 1) + (x_max1 - day_prop_aligned.index(t_dates[len(t_dates) - 1]))
                else:
                    x_max = list(t_dates).index(str(sys.argv[6])) 
            except:
                x_min = None
        
            if (x_min is not None):
                t_dates_show = np.array(t_dates)[x_min:x_max+1]
                check_last = x_max1
            else:
                t_dates_show = t_dates
                check_last = len(day_prop_aligned) - 1
                
            if len(t_dates_show)>200:
                pp = 7*4
            else:
                pp = min(len(t_dates_show), 14)
            
            perday = np.arange(0,len(t_dates_show), pp)
            date_ticks = t_dates_show[perday].tolist()
            if t_dates[len(t_dates) - 1] not in date_ticks:
                n=list(t_dates).index(date_ticks[-1])+pp
                while n<len(t_dates)-1:
                    date_ticks.append(t_dates[n])
                    perday = np.append(perday, n)
                    n += pp
                date_ticks.append(t_dates[len(t_dates) - 1])
                perday = np.append(perday, len(t_dates) - 1)
            
            change = len(date_ticks)
            
            if day_prop_aligned[check_last] not in date_ticks:
                try:
                    n=list(day_prop_aligned).index(date_ticks[-1])+pp
                except:
                    n= perday[-1]+pp
                while n<len(day_prop_aligned)-1:
                    date_ticks.append(day_prop_aligned[n])
                    perday = np.append(perday, n)
                    n += pp
                if check_last-perday[-1]<np.ceil(pp/5):
                    perday = perday[:-1]
                    date_ticks = date_ticks[:-1]
                date_ticks.append(day_prop_aligned[check_last])
                perday = np.append(perday, check_last)
              
            perday_orig = []
            for i in range(len(np.array(date_ticks)[:change])):
                try:
                    perday_orig.append(list(t_dates).index(date_ticks[i]))
                except:
                    perday_orig.append(perday[i])
            
            for j in range(len(np.array(date_ticks[change:]))):
                try:
                    perday_orig.append(list(day_prop_aligned).index(date_ticks[change+j]))
                except:
                    perday_orig.append(perday[change+j])
            
            ax_k.set_xticks(perday_orig)
            ax_k.set_xticklabels(date_ticks,
                rotation = 45, horizontalalignment = "right")
            
            if (x_min is not None):
                ax_k.set_xlim((x_min, x_max))
                ax_k_twin.set_xlim((x_min1, x_max1))
                
            #ax_twin.set_ylim((-0.02, 0.02))
            ax_k.axhline(xmin = 0, xmax = len(t_prop_aligned), ls = "--", linewidth = 2, color = "black")
            ax_k.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
            ax_k_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
            ax_k.set_ylabel("Relative fitness", fontsize = 20)
            ax_k_twin.set_ylabel("Spikegroups Frequency (daily %)", fontsize = 20)
            pdf_k = PdfPages(sys.argv[w_save]+"/relative_fitness_%s_vs_prop.pdf"%(lab_k_fn))
            pdf_k.savefig(fig_k, bbox_inches = "tight")
            pdf_k.close()
     
            fig_k.savefig(sys.argv[w_save]+"/relative_fitness_%s_vs_prop.svg"%(lab_k_fn), bbox_inches = "tight")
            plt.close()
        
            ### Separate figure for relative fitness vs change in proportion
            PreFig(xsize = 20, ysize = 20)
            fig2 = plt.figure(figsize = (15, 7))
            ax2 = fig2.add_subplot(1, 1, 1)
            # different axis for proportions
            ax2_twin = ax2.twinx()
            
            gamma_prop_masked = ma.masked_array(gamma_prop, mask = SI_mask)
            ax2.fill_between(inds_dates, gamma_SI_min, gamma_SI_max, color = "green", alpha = 0.3, label = lab_k)
            ax2_twin.plot(inds_dates, gamma_prop_masked, color = "orange", label=lab_k)
            #ax2_twin.scatter(inds_dates, gamma_prop_masked, marker = ".", color = "orange")

            ax2.set_xticks(perday_orig)
            ax2.set_xticklabels(date_ticks,
                rotation = 45, horizontalalignment = "right")
            
            if x_min is not None:
                ax2.set_xlim((x_min, min(x_max, inds_dates[-1])))
                ax2_twin.set_xlim((x_min, min(x_max, inds_dates[-1])))   
                
            
            ymin1, ymax1 = ax2.get_ylim()
            ymin2, ymax2 = ax2_twin.get_ylim()
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            
            #loc0 = min(np.abs(ymin1)/(np.abs(ymin1)+np.abs(ymax1)), np.abs(ymax1)/(np.abs(ymin1)+np.abs(ymax1)))
            #mpl_axes_aligner.align.yaxes(ax2, 0, ax2_twin, 0, loc0)
            mpl_axes_aligner.align.yaxes(ax2, 0, ax2_twin, 0, 0.5)
            
            if (ymin1/ymin2 >0.5) or (ymax1/ymax2>0.5) or (ymin2/ymin1 >0.5) or (ymax2/ymax1>0.5):
                ax2.set_ylim((ymin, ymax))
                ax2_twin.set_ylim((ymin, ymax))   
                
            ax2.axhline(xmin = 0, xmax = len(day_prop_aligned), ls = "--", linewidth = 2, color = "black")
            ax2.set_ylabel("Relative fitness $\gamma_y$", fontsize = 20)
            ax2_twin.set_ylabel("Change in proportion $\gamma_{prop}$", fontsize = 20)
            ax2.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
            ax2_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
            pdf_2 = PdfPages(sys.argv[w_save]+"/relative_fitness_%s.pdf"%lab_k_fn)
            pdf_2.savefig(fig2, bbox_inches = "tight")
            pdf_2.close()
     
            fig2.savefig(sys.argv[w_save]+"/relative_fitness_%s.svg"%lab_k_fn, bbox_inches = "tight")
            plt.close()
            
            status_list.append(lab_status)

        else:
            print("No lineage in group %s has E[Susceptible] available, if needed, first compute it in main config"%lineage_list[k])
            status_list.append("No data")
    
    
    ax.axhline(xmin = 0, xmax = len(day_prop_aligned), ls = "--", linewidth = 2, color = "black")
    
    ymin1, ymax1 = ax.get_ylim()
    ymin2, ymax2 = ax_twin.get_ylim()
    #ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
    ax.set_ylim((ymin1, ymax1))
    #ax_twin.set_ylim((ymin2, ymax2))
    #loc0 = min(np.abs(ymin1)/(np.abs(ymin1)+np.abs(ymax1)), np.abs(ymax1)/(np.abs(ymin1)+np.abs(ymax1)))
    #mpl_axes_aligner.align.yaxes(ax, 0, ax_twin, 0, loc0)
    try:
        if str(sys.argv[5]) in list(t_dates):
            x_min = list(t_dates).index(str(sys.argv[5]))
        else:
            x_min = 0
        if str(sys.argv[5]) in list(day_prop_aligned):
            x_min1 = day_prop_aligned.index(str(sys.argv[5]))
        else:
            x_min1 = 0
        if str(sys.argv[6]) in list(day_prop_aligned):
            x_max1 = day_prop_aligned.index(str(sys.argv[6]))
        else:
            x_max1 = len(day_prop_aligned) - 1
            
        if str(sys.argv[6]) not in list(t_dates):
            x_max = (len(t_dates) - 1) + (x_max1 - day_prop_aligned.index(t_dates[len(t_dates) - 1]))
        else:
            x_max = list(t_dates).index(str(sys.argv[6]))
    except:
        x_min = None

    if (x_min is not None):
        t_dates_show = np.array(t_dates)[x_min:x_max+1]
        check_last = x_max1
    else:
        t_dates_show = t_dates
        check_last = len(day_prop_aligned) - 1

    if len(t_dates_show)>200:
        pp = 7*4
    else:
        pp = min(len(t_dates_show), 14)
    
    perday = np.arange(0,len(t_dates_show), pp)
    date_ticks = t_dates_show[perday].tolist()
    if t_dates[len(t_dates) - 1] not in date_ticks:
        n=list(t_dates).index(date_ticks[-1])+pp
        while n<len(t_dates)-1:
            date_ticks.append(t_dates[n])
            perday = np.append(perday, n)
            n += pp
        date_ticks.append(t_dates[len(t_dates) - 1])
        perday = np.append(perday, len(t_dates) - 1)
    
    change = len(date_ticks)
    
    if day_prop_aligned[check_last] not in date_ticks:
        try:
            n=list(day_prop_aligned).index(date_ticks[-1])+pp
        except:
            n= perday[-1]+pp
            
        while n<len(day_prop_aligned)-1:
            date_ticks.append(day_prop_aligned[n])
            perday = np.append(perday, n)
            n += pp
            
        if check_last-perday[-1]<np.ceil(pp/5):
            perday = perday[:-1]
            date_ticks = date_ticks[:-1]
        date_ticks.append(day_prop_aligned[check_last])
        perday = np.append(perday, check_last)
    
    perday_orig = []
    for i in range(len(np.array(date_ticks)[:change])):
        try:
            perday_orig.append(list(t_dates).index(date_ticks[i]))
        except:
            perday_orig.append(perday[i])
   
    for j in range(len(np.array(date_ticks[change:]))):
        try:
            perday_orig.append(list(day_prop_aligned).index(date_ticks[change+j]))
        except:
            perday_orig.append(perday[change+j])
    
    ax.set_xticks(perday_orig)
    ax.set_xticklabels(date_ticks,
        rotation = 45, horizontalalignment = "right")
    
    if (x_min is not None):
        ax.set_xlim((x_min, x_max))
        ax_twin.set_xlim((x_min1, x_max1))
    
    #ax.legend(loc = (1.2, 0.) ,fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.27),fancybox=True, shadow=True, ncol=5)
    #ax_twin.legend(loc = (1.2, 0.), fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
    ax.set_ylabel("Relative fitness", fontsize = 20)
    ax_twin.set_ylabel("Spikegroup Frequency (monthly)", fontsize = 20)
    pdf = PdfPages(sys.argv[w_save]+"/relative_fitness_groups.pdf")
    pdf.savefig(fig, bbox_inches = "tight")
    pdf.close()
    fig.savefig(sys.argv[w_save]+"/relative_fitness_groups.svg", bbox_inches = "tight")
    plt.close()

    
    ### work on spikes group props fig
    try:
        if str(sys.argv[5]) in list(day_prop_aligned):
            x_min1 = day_prop_aligned.index(str(sys.argv[5]))
        else:
            x_min1 = 0
        if str(sys.argv[6]) in list(day_prop_aligned):
            x_max1 = day_prop_aligned.index(str(sys.argv[6]))
        else:
            x_max1 = len(day_prop_aligned) - 1
    except:
        x_min1 = None
    
    if (x_min1 is not None):
        t_show = np.array(day_prop_aligned)[x_min1:x_max1+1]
        check_last = x_max1
    else:
        t_show = day_prop_aligned
        check_last = len(day_prop_aligned) - 1
    
    perday = np.arange(0,len(t_show), pp)
    date_ticks = np.array(t_show)[perday].tolist()
    
    if day_prop_aligned[check_last] not in date_ticks:
        try:
            n=list(day_prop_aligned).index(date_ticks[-1])+pp
        except:
            n=perday[-1] + pp
            
        while n<len(day_prop_aligned)-1:
            date_ticks.append(day_prop_aligned[n])
            perday = np.append(perday, n)
            n += pp
        if check_last-perday[-1]<np.ceil(pp/5):
            perday = perday[:-1]
            date_ticks = date_ticks[:-1]
        
        date_ticks.append(day_prop_aligned[check_last])
        perday = np.append(perday, check_last)
       
    
    perday_orig = []
    for i in range(len(np.array(date_ticks))):
        try:
            perday_orig.append(list(day_prop_aligned).index(date_ticks[i]))
        except:
            perday_orig.append(perday[i])
            
    ax_prop.set_xticks(perday_orig)
    ax_prop.set_xticklabels(date_ticks,
        rotation = 45, horizontalalignment = "right")
    
    if (x_min1 is not None):
        ax_prop.set_xlim((x_min1, x_max1))
    
    already_prop = ma.masked_array(already_prop_aligned, mask=prop_mask_aligned)
    return status_list, already_prop, ax_prop, perday_orig, fig_prop, t_prop_aligned

num_groups = int(sys.argv[7])
w_save = 8
k = 9
lineage_list = []
color_list = []
custom_col = sns.color_palette("Set2", 100) 
s = 0

for i in range(num_groups):
    lineage_list.append(str(sys.argv[k+i]))
    try:
        if "/" not in str(sys.argv[k+num_groups+i]):
            color_list.append(str(sys.argv[k+num_groups+i]))
        else:
            split_col = str(sys.argv[k+num_groups+i]).split("/")
            color_list.append(tuple([float(split_col[c]) for c in range(len(split_col))])) ### anything else is error)
    except:
        rand_num = np.random.choice(1, 100)
        if s<len(custom_col):
            color_list.append(custom_col[s])
        else:
            color_list.append(sns.color_palette("rocked", rand_num)[0])
        s +=1

status_list, already_prop, ax_prop, perday_orig, fig_prop, t_prop_aligned = plot_fit(ES_lin_dir, lineage_list, color_list, w_save, already_prop = np.zeros((len(t_prop))))
### Group Plot proportion of all other spikegroups
ax_prop.plot(t_prop_aligned, (100 - 100*already_prop), linewidth = 3, color = "grey", label = "Other")
#ax_prop.scatter(t_prop_aligned, (100 - 100*already_prop), marker = ".", color = "grey")
ymin, ymax = ax_prop.get_ylim()

ax_prop.set_ylim(((0, 1.0*ymax)))
ax_prop.legend(loc = (1.2, 0.), fontsize = 20, ncols = np.ceil(len(lineage_list)/4).astype(int))
pdf2 = PdfPages(sys.argv[w_save]+"/Groups_proportions.pdf")
ax_prop.set_ylabel("Spikegroups Frequency (daily %)", fontsize = 20)
pdf2.savefig(fig_prop, bbox_inches = "tight")
fig_prop.savefig(sys.argv[w_save]+"/Groups_proportions.svg")
pdf2.close()

status = pd.DataFrame({"lineage":lineage_list, "spikegroups_found":status_list})
status.to_csv(sys.argv[w_save]+"/plot_status.csv")
 

        