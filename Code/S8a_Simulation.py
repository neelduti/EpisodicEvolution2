#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages


def MakeOrChangeDir(Dr, change=True):
    if os.path.isdir(Dr): pass
    else: os.mkdir(Dr)
    if change: os.chdir(Dr)


EventDf = pd.read_csv('../CommonDf.tsv', sep='\t', 
                       usecols=['HS', 'HS_Subs', 'PT_Subs', '#1#_Subs',
                                'GG_Subs', 'NL_Subs',], 
                       index_col = 'HS')

#change the current working directory and subdirectories
MakeOrChangeDir('Simulation')
MakeOrChangeDir('Suppli', change=False)


EventDf['Events'] = EventDf[['HS_Subs', 'PT_Subs', 
                                'GG_Subs', 'NL_Subs', '#1#_Subs',]].sum(axis=1)
EventDf = EventDf[EventDf['Events'] > 0]
EventDf = EventDf.astype(int)
a = ['h', 'c', 'p', 'g', 'n']

PvalDf = pd.DataFrame(columns=a)

n=0
for Events in EventDf.Events.value_counts().sort_index().index:

    OutcomeDf = pd.DataFrame(columns=a)
    for k in range(100000):
        out = pd.Series(np.random.choice(a, Events, p=[0.1174, 0.1115,  0.0372, 0.1542, 0.5797])).value_counts()
        OutcomeDf.loc[k, out.index] = out.values
    OutcomeDf.replace(np.nan, 0, inplace=True)
    
    ValueCountDf = pd.DataFrame(index = range(Events + 1))
    for x in a: 
        ValueCountDf[x] = OutcomeDf[x].value_counts()
    ValueCountDf.replace(np.nan, 0, inplace=True)
    ValueCountDf.rename(columns={'h':'Human', 'c':'Chimpanzee', 'p':'#1#', 'g':'Gorilla', 'n':'Gibbon'}, inplace=True)
    
    ValueCountDf.to_csv(f'Suppli/ValueCountDf_{Events:0>3d}_Events.tsv', sep='\t', header=True, index=True, index_label='Events')
    
    ValueCountDf['Events'] = ValueCountDf.index
    with plt.rc_context(dict(sns.axes_style("white", rc = {
    'axes.spines.right': False,
    'axes.spines.top': False,}),
                         **sns.plotting_context("paper", rc={
                             'font.size': 14.0,
                             'axes.labelsize':14, 
                             'axes.titlesize': 20,
                             'xtick.labelsize': 14,
                             'ytick.labelsize': 14,
                             'legend.fontsize': 12,
                             }))):    
        plt.close('all')
        plt.clf()
        plt.rcParams['pdf.fonttype'] = 'truetype'
        fig2 = plt.figure(figsize=(10,5), dpi=400)
        spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig2)
        ax1 = fig2.add_subplot(spec2[0, 0])
    
        sns.lineplot(x='Events', y="Frequency", hue='Branches', 
                     data = ValueCountDf.melt('Events', var_name = 'Branches', value_name = "Frequency"), 
                     ax=ax1)
        
        ax1.set_title(f'Events = {Events}', loc='Center', fontweight='bold') 
        plt.tight_layout(pad=1.08)
        plt.savefig(f'Suppli/Events_{Events:0>3d}_distribution.png')
       
    #create p-value confidence interval through two sided rank test ranking
    for Indx,row in EventDf[EventDf.Events == Events].iterrows(): 
        Observed = pd.Series({'h':row['HS_Subs'], 'c':row['PT_Subs'], 'p':row['#1#_Subs'], 'g':row['GG_Subs'], 'n':row['NL_Subs']})
        OutcomeDf.loc['obs', Observed.index] = Observed.values
        RankDf = pd.DataFrame(index=a)
        RankDf['rank'] = OutcomeDf.rank().loc['obs',:]-1
        RankDf['Emp1'] = ((RankDf['rank'])/(len(OutcomeDf)-1))
        RankDf['Emp2'] = 1 - RankDf['Emp1']
        RankDf['P_value'] = RankDf[['Emp1', 'Emp2']].min(axis=1)*2
        PvalDf.loc[Indx, RankDf.P_value.index] = RankDf.P_value.values
    
    # n+=1
    # if n > 20: break

PvalDf.rename(columns={'h':'HS', 'c':'PT', 'p':'#1#', 'g':'GG', 'n':'NL'}, inplace=True)
PvalDf.loc[PvalDf.index, 'Events'] = EventDf.loc[PvalDf.index, 'Events']
PvalDf.to_csv('RankTestPvalue.tsv', sep='\t', header=True, index=True, index_label='Id')
