import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import pdb
import mpl_axes_aligner

def moving_average(X, window = 7):
    u = np.zeros(len(X))
    u[:window] = X[:window]
    for i in range(window, len(X)):
        u[i] = np.mean(X[i-window:i+1])
    
    return u

lineage_freq = pd.read_csv(sys.argv[2])
threshold = float(sys.argv[3])
variant = str(sys.argv[4])
S_mean_file = sys.argv[5]


#ES_df = pd.read_csv("demo/results/Immunological_Landscape/Susceptible_SpikeGroup_lineage_XXX_all_PK.csv")
#lineage_freq = pd.read_csv("demo/results/Daily_Lineages_Freq.csv")
#threshold = 0 #(percent)
#variant = "BA.5.1"

# needs to be updated to allow individual weighting 
S_mean_df = pd.read_csv(S_mean_file)
S_all_mean = S_mean_df.to_numpy()[:, (S_mean_df.columns != "Days")&(S_mean_df.columns != "Unnamed: 0")].astype(float)
t_dates = S_mean_df["Days"].tolist()

# processing of frequency data
try:
    lineage_freq.drop(columns = "Unnamed: 0", inplace = True)
except:
    pass

lineage_freq = lineage_freq[lineage_freq['date'].isin(t_dates)]
prop_mask = np.all(lineage_freq.loc[:, lineage_freq.columns != 'date'] == 0, axis = 1)

freqs = lineage_freq.loc[:, lineage_freq.columns != 'date']
# imputing frequencies below threshold and normalization
freqs = freqs.mask(freqs < threshold)
freqs = freqs.fillna(0)
col_sums = freqs.sum(axis = 1).values
freqs = freqs.divide(col_sums, axis="rows")
freqs = freqs.fillna(0)
lineage_freq.loc[:, lineage_freq.columns != 'date'] = freqs
day_prop = lineage_freq["date"].tolist()

import matplotlib
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)
    
def plot_fit(ES_df, lineage, w_save = 6):
    # processing of susceptibles
    try:
        ES_df.drop(columns = "Unnamed: 0", inplace = True)
    except:
        pass
    
    es_cols = ES_df.columns
    ES_df = ES_df[ES_df['Days'].isin(t_dates)]
    ES_ranges = ES_df.to_numpy()[:, es_cols!="Days"].astype(float)
    
    # calculation of change in relative frequency from model
    gamma_SI = np.zeros((len(t_dates), ES_ranges.shape[1]))
    
    for i in range(ES_ranges.shape[1]):
        S_x = ES_ranges[:, i]
        
        S_mean = S_all_mean[:, i]
        gamma_SI[:, i] = np.divide(S_x - S_mean, S_mean, out = S_x, where = S_mean != 0)
        
    # get min max gamma over PK at each timepoints
    gamma_SI_min, gamma_SI_max = np.min(gamma_SI, axis = 1), np.max(gamma_SI, axis = 1)
    
    # change in relative frequency from genomic surveillance data 
    if "Spike. " + lineage in lineage_freq.columns.astype(str):
        Pseudo_Prop = moving_average(lineage_freq["Spike. " + lineage], window = 14)
        Pseudo_Prop[Pseudo_Prop < 0.05] = 0        
        #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
    elif lineage in lineage_freq.columns.astype(str):
        Pseudo_Prop = moving_average(lineage_freq[lineage], window = 14)
        Pseudo_Prop[Pseudo_Prop < 0.05] = 0
        #Pseudo_Prop = Pseudo_Prop/np.sum(Pseudo_Prop)
    else:
        Pseudo_Prop = np.zeros(len(t_dates))
    
    gamma_prop = np.zeros(len(t_dates))
    Pseudo_Prop = list(Pseudo_Prop)
    SI_mask = np.zeros(len(t_dates)).astype(bool)
    for l in range(len(t_dates)):
        if t_dates[l] in day_prop:
            w_l = list(day_prop).index(t_dates[l])
            if w_l + 1 < len(Pseudo_Prop):
                if Pseudo_Prop[w_l] == 0 or Pseudo_Prop[w_l+1] == 0:
                    gamma_prop[l] = float('nan')
                    #SI_mask[l] = True
                else:
                    gamma_prop[l] = (Pseudo_Prop[w_l+1]/Pseudo_Prop[w_l]) - 1
            else:
                gamma_prop[l] = float("nan")
                SI_mask[l] = True
        else:
            gamma_prop[l] = float('nan')
            SI_mask[l] = True
    
    # plotting
    PreFig(xsize = 20, ysize = 20)
    fig = plt.figure(figsize = (15, 7))
    ax = fig.add_subplot(1, 1, 1)
    ax_twin = ax.twinx()
    
    SI_mask = np.array(SI_mask) + prop_mask[:len(t_dates)]
    gamma_SI_min = ma.masked_array(gamma_SI_min, mask = SI_mask)
    gamma_SI_max = ma.masked_array(gamma_SI_max, mask = SI_mask)
    gamma_prop = ma.masked_array(gamma_prop, mask = SI_mask)
    ax.fill_between(t_dates, gamma_SI_min, gamma_SI_max, color = "green", alpha = 0.3, label = "$\gamma_{%s}$"%lineage)
    #more smoothing
    #gamma_prop = moving_average(gamma_prop, window = 14)
    ax_twin.plot(t_dates, gamma_prop, color = "orange", linewidth = 4, label="$\gamma_prop$ %s"%lineage)
    
    #ax.axhline(xmin = 0, xmax = t_dates[-1], ls = "--", linewidth = 2, color = "black")
    ax.axhline(xmin = 0, xmax = len(t_dates), ls = "--", linewidth = 2, color = "black")
    
    try:
        if str(sys.argv[7]) in list(t_dates):
            x_min = list(t_dates).index(str(sys.argv[7]))
        else:
            x_min = 0
        if str(sys.argv[8]) not in list(t_dates):
            x_max = (len(t_dates) - 1) 
        else:
            x_max = list(t_dates).index(str(sys.argv[8])) 
    except:
        x_min = None
        x_max = None
    
    if (x_min is not None):
        ax.set_xlim((x_min, x_max))
        t_dates_show = np.array(t_dates)[x_min:x_max+1]
    else:
        t_dates_show = t_dates

    if len(t_dates_show)>200:
        pp = 7*4
    else:
        pp = min(len(t_dates_show), 14)

    perday = np.arange(0,len(t_dates_show), pp)
    date_ticks = t_dates_show[perday].tolist()
    if t_dates[len(t_dates) - 1] not in date_ticks:
        try:
            n=list(t_dates).index(date_ticks[-1])+pp
            while n<len(t_dates)-1:
                date_ticks.append(t_dates[n])
                perday = np.append(perday, n)
                n += pp
            date_ticks.append(t_dates[len(t_dates) - 1])
            perday = np.append(perday, len(t_dates) - 1)
        except:
            pass
        
    if x_min is not None:
        perday_orig = []
        for i in range(len(date_ticks)):
            perday_orig.append(list(t_dates).index(date_ticks[i]))
    else:
        perday_orig = perday
        
    ax.set_xticks(perday_orig)
    ax.set_xticklabels(date_ticks,
        rotation = 45, horizontalalignment = "right")
    
    if lineage == "EG.1.3":
        print(4)
        pdb.set_trace() 
    
    ymin1, ymax1 = ax.get_ylim()
    ymin2, ymax2 = ax_twin.get_ylim()
    ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
    mpl_axes_aligner.align.yaxes(ax, 0, ax_twin, 0, 0.5)
    # Align y = 0 of ax1 and ax2 with the center of figure.
    #loc0 = min(np.abs(ymin)/(np.abs(ymin)+np.abs(ymax)), np.abs(ymax)/(np.abs(ymin)+np.abs(ymax)))
    #mpl_axes_aligner.align.yaxes(ax, 0, ax_twin, 0, loc0)
    if (ymin1/ymin2 >0.5) or (ymax1/ymax2>0.5) or (ymin2/ymin1 >0.5) or (ymax2/ymax1>0.5):
        ax.set_ylim((ymin, ymax))
        ax_twin.set_ylim((ymin, ymax))        
    
    ax.set_ylabel("Relative fitness $\gamma_y$", fontsize = 20)
    ax_twin.set_ylabel("Change in proportion $\gamma_{prop}$", fontsize = 20)
    
    pdf = PdfPages(sys.argv[w_save]+"/relative_fitness_%s.pdf"%lineage)
    pdf.savefig(fig, bbox_inches = "tight")
    pdf.close()
 
    fig.savefig(sys.argv[w_save]+"/relative_fitness_%s.svg"%lineage, bbox_inches = "tight")
    plt.close()
    
if variant != "ALL":
    try:
        ES_df = pd.read_csv(sys.argv[1])
    except:
        try:
            ES_df = pd.read_csv(sys.argv[1]+"/Susceptible_SpikeGroup_%s_all_PK.csv"%variant)
        except:
            print("Computation needed: Expected Susceptible file is not available for %s"%variant)
    w_save = 6
    plot_fit(ES_df, variant, w_save = w_save)
    status = pd.DataFrame({"lineage":variant, "relative_advantage":"Done"}, index = [1])
    status.to_csv(sys.argv[6]+"/plot_status.csv")

else:
    status_list = []
    lineage_freq.drop(columns = "date", inplace = True)
    w_save = 6
    num = 1
    for i in range(len(lineage_freq.columns.astype(str))):
        variant = lineage_freq.columns.astype(str)[i]
        print("Plot relative fitness of %s (%d/%d)"%(variant, i+1, len(lineage_freq.columns.astype(str))))
        try:
            if "Spike" in variant:
                ES_df = pd.read_csv(sys.argv[1]+"/Susceptible_SpikeGroup_%s_all_PK.csv"%variant[7:])
            else:
                ES_df = pd.read_csv(sys.argv[1]+"/Susceptible_SpikeGroup_%s_all_PK.csv"%variant)
            
            plot_fit(ES_df, variant, w_save = w_save)
            status_list.append("Done")
        except:
            print(num, "Was not computed: %s Not present in file Cross_react_dic_spikegroups_ALL.pck"%variant[7:])
            num +=1
            status_list.append("Not computed: %s Not present in file Cross_react_dic_spikegroups_ALL.pck"%variant[7:])

    status = pd.DataFrame({"lineage":lineage_freq.columns.astype(str), "relative_advantage":status_list})
    status.to_csv(sys.argv[w_save]+"/plot_status_all.csv")
        

