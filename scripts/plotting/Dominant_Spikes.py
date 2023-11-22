#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 15:39:36 2023

@author: raharinirina
"""

import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

lineage_freq = pd.read_csv(sys.argv[1])

try:
    lineage_freq = lineage_freq.drop("Unnamed: 0", axis = 1)
except:
    pass

try:
    lineage_freq = lineage_freq.drop("date", axis = 1)
except:
    pass

try:
    lineage_freq = lineage_freq.drop("week_num", axis = 1)
except:
    pass


import matplotlib
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)

PreFig(xsize = 20, ysize = 20)
fig = plt.figure(figsize = (9, 7))
ax = fig.add_subplot(1, 1, 1)

### Sort lineage_freqs
sorted_cols = lineage_freq.columns[np.argsort(np.sum(lineage_freq.to_numpy().astype(float),  axis = 0))[::-1]] ## sort in descending order of total proportions over the entire timeframe
cols = sns.color_palette("husl", len(sorted_cols)) 
count = 0
for spike in sorted_cols:
    if max(lineage_freq[spike]) > 0:
        ax.plot(lineage_freq[spike], label = spike, linewidth = 3, color = cols[count])
        count +=1

try:
    ax.legend(loc = (1.2, 0), ncols = (np.ceil(count)//15).astype(int))
except:
    ax.legend(loc = (1.2, 0), ncols = (np.ceil(count)//15).astype(int))

ax.set_title("Ordered Spikesgroups by total sum of proportions", fontsize = 20)
pdf2 = PdfPages(sys.argv[2]+"/SpikeGroups_Props_overview.pdf")
ax.set_ylabel(sys.argv[2]+"Frequency (daily %)", fontsize = 20)
pdf2.savefig(fig, bbox_inches = "tight")
pdf2.close()
fig.savefig(sys.argv[2]+"/SpikeGroups_Props_overview.svg")
