#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
import numpy as np
import os

pd.set_option('display.max_columns', 15)
pd.set_option('display.width', 200)
sns.set(context='paper')
sns.set_palette("colorblind")

#Load the original orthlog table
CommonDf = pd.read_csv('../CommonDf.tsv', sep='\t', index_col=0) #12621
#drop gene families with no overlap
CommonDf.dropna(inplace=True, subset=['overlap']) #12618
AaLs = [x for x in list(CommonDf.columns) if x.endswith('_aa')]


CommonDf['Min'] = CommonDf[AaLs].apply(lambda row: min(row), axis=1)
CommonDf['Max'] = CommonDf[AaLs].apply(lambda row: max(row), axis=1)

################# figure 2
from scipy.stats import linregress
linregress(CommonDf.Min, CommonDf.Max)

CommonDf['MaxDif'] = CommonDf[AaLs].apply(lambda row: (max(row) - min(row)), axis=1)
CommonDf['RelMaxDif'] = CommonDf['MaxDif']*100/CommonDf['Min']

CommonDf['RelOverlap'] = CommonDf['overlap']*100/CommonDf['Min']
CommonDf['%MinId'] = CommonDf['AbsId']*100/CommonDf['Min']

with plt.rc_context(dict(sns.axes_style("white", rc = {
    'axes.spines.right': False,
    'axes.spines.top': False,}),
                         **sns.plotting_context("paper", rc={
                             'axes.labelsize':14, 
                             'axes.titlesize': 18,
                             'xtick.labelsize': 12,
                             'ytick.labelsize': 12,}))):
    plt.close('all')
    plt.clf()
    plt.rcParams['pdf.fonttype'] = 'truetype'
    fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(nrows=1, ncols=4,figsize=(10,4))
    #figure 2a scatterplot for min and max distribution of each cluster
    sns.regplot(x="Min", y="Max", data=CommonDf, y_jitter=0.05, line_kws={'color':'red'}, ax=ax1)
    ax1.set(xticks=range(0,8001,2000), xlim=(0,8000), ylim=(0,8000), xlabel="Shortest Ortholog\n(aa residues)", ylabel="Longest Ortholog\n(aa residues)")
    plt.setp(ax1.get_xmajorticklabels(), rotation=30)
    ax1.text(-0.3, 1.2, 'a', transform=ax1.transAxes,
      fontsize=30, fontweight='bold', va='top', ha='right')
    #figure 2b %maxdiff histogram
    sns.histplot(CommonDf[CommonDf['RelMaxDif'] > 0]['RelMaxDif'], bins=range(0,15,1), kde=False, ax=ax2)
    #counts, bins, patches = ax3.hist(CommonDf[CommonDf['%MaxDif'] > 0]['%MaxDif'], bins=10)
    ax2.set(xticks=range(1,15,1), xlim=(0,9), xlabel="Length difference (%)", ylabel="Ortholog Families")
    plt.setp(ax2.get_xmajorticklabels(), rotation=30)
    ax2.text(-0.3, 1.2, 'b', transform=ax2.transAxes,
      fontsize=30, fontweight='bold', va='top', ha='right') 
    # figure 2c Rel-overlap
    sns.histplot(CommonDf['RelOverlap'], bins=range(90,101,1), kde=False, ax=ax3)
    ax3.set(xticks=range(90,101,2), xlim=(90,100), xlabel="Aligment Saturation (%)", ylabel="Ortholog Families")
    plt.setp(ax3.get_xmajorticklabels(), rotation=30)
    ax3.text(-0.3, 1.2, 'c', transform=ax3.transAxes,
      fontsize=30, fontweight='bold', va='top', ha='right') 
    # figure 2d %Abs identiy
    sns.histplot(CommonDf['%AbsId'], bins=range(90,101,1), kde=False, ax=ax4)
    ax4.set(xticks=range(90,101,2), xlim=(90,100), xlabel="Idenity (%)\n(within aligned region)", ylabel="Ortholog Families")
    plt.setp(ax4.get_xmajorticklabels(), rotation=30)
    ax4.text(-0.3, 1.2, 'd', transform=ax4.transAxes,
      fontsize=30, fontweight='bold', va='top', ha='right') 
    plt.tight_layout(pad=1.08)
    plt.savefig('Fig2_OrthologFamilies.pdf')
    plt.savefig('Fig2_OrthologFamilies.png')

#Min in sp
for Index, Row in CommonDf.iterrows():
    m = 0
    for Sp in AaLs:
        if Row[Sp] == Row['Min']:
            CommonDf.at[Index, 'Min_sp'] = Sp
            m += 1
    if m > 1: CommonDf.at[Index, 'Min_sp'] = m
CommonDf.Min_sp.value_counts()

for Index, Row in CommonDf.iterrows():
    m = 0
    for Sp in AaLs:
        if Row[Sp] == Row['Max']:
            CommonDf.at[Index, 'Max_sp'] = Sp
            m += 1
    if m > 1: CommonDf.at[Index, 'Max_sp'] = m
CommonDf.Max_sp.value_counts()

################# figure 3 Multi sp substitution
CommonDf['MultiSp_Subs'] = (CommonDf['Subs'] - CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs']].sum(axis=1))
CommonDf['%MultiSp_Subs'] = (CommonDf['Subs'] - CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs']].sum(axis=1))*100/CommonDf['overlap']

plt.close('all')
plt.clf()
plt.rcParams['pdf.fonttype'] = 'truetype'
fig, ax2 = plt.subplots(nrows=1, ncols=1,figsize=(3,5))
# figure 3g sum of all substitutions
CommonDf[['#1#_Subs', 'Convergent_Subs', 'OneInOutId', 'OnlyInGpId', 'OnlyOutGpId',
       'NoId',]].sum().plot.bar(ax=ax2)
ax2.set(ylabel="Total Number of sites", xticklabels=["#1#_Subs", 'Two identities - inconsistent', 'One identity - inconsistent', 'Only Ingroup identical', 'Only Outgroup identical',
       'No identity',])
plt.setp(ax2.get_xmajorticklabels(), rotation=30, ha='right')

plt.tight_layout(pad=1.08)
plt.savefig('Fig3g_MultiSpSubs.pdf')

################# figure 4 Departure from the average tree
with plt.rc_context(dict(sns.axes_style("white", rc = {
    'axes.spines.right': False,
    'axes.spines.top': False,}),
                         **sns.plotting_context("paper", rc={
                             'font.size': 14.0,
                             'axes.labelsize':14, 
                             'axes.titlesize': 40,
                             'xtick.labelsize': 14,
                             'ytick.labelsize': 14,
                             'legend.fontsize': 12,
                             }))):    
    plt.close('all')
    plt.clf()
    plt.rcParams['pdf.fonttype'] = 'truetype'
    fig2 = plt.figure(figsize=(10,5))
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig2)
    #ax1 = fig2.add_subplot(spec2[0, 0])
    ax2 = fig2.add_subplot(spec2[0, 0])
    #1 RF metric
    PcDf = pd.read_csv('../Subs_gblocks.aln/WithNull/PoissonCorrectedDf.tsv', sep='\t', index_col=0)
    BraScStatSigDf = pd.read_csv('../Subs_gblocks.aln/WithNull/PoissonCorrectedDf/FdrBraScStatSigDf.tsv', sep='\t', index_col=0)
    BraScDf = pd.read_csv('../Subs_gblocks.aln/WithNull/PoissonCorrectedDf/BranchScoreDf.tsv', sep='\t', index_col=0)
    Insig = BraScDf.loc[BraScStatSigDf[BraScStatSigDf.TreeFull >= 0.05].index, :].TreeFull
    Sig = BraScDf.loc[BraScStatSigDf[BraScStatSigDf.TreeFull < 0.05].index, :].TreeFull

    #2 Histogram Tree length
    iSign = PcDf.loc[BraScStatSigDf[BraScStatSigDf.TreeFull < 0.05].index.intersection(PcDf.index), ['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']].sum(axis=1).sort_values(ascending=False)
    iInSign = PcDf.loc[BraScStatSigDf[BraScStatSigDf.TreeFull >= 0.05].index.intersection(PcDf.index), ['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']].sum(axis=1).sort_values(ascending=False)
    ax2.hist([iInSign, iSign], bins=np.linspace(0,0.07,71), stacked=True, label=['Insignificant', 'Significant'], color=['grey', 'black'])
    MyMean = PcDf.loc[BraScStatSigDf.index.intersection(PcDf.index), ['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']].sum(axis=1).mean()
    ax2.axvline(MyMean, color='k',linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = ax2.get_ylim()
    ax2.text(MyMean*1.05, max_ylim*0.9, 'Mean:\n{:.3f}'.format(MyMean))
    ax2.set(xlabel='PC Tree Length', ylabel = 'Ortholog Families')
    ax2.legend(title='Deviation from the\nmean tree', loc='upper right', title_fontsize=12, fancybox=True, framealpha=0.5)
    # ax2.set_title('a', loc='left', fontweight='bold') 
    # ax2.text(-0.15, 1.1, 'b', transform=ax2.transAxes,
    #   fontsize=40, fontweight='bold', va='top', ha='right')
    plt.tight_layout(pad=1.08)
    plt.savefig('Fig4_TreeDeparture.pdf')
    plt.savefig('Fig4_TreeDeparture.png')

for Cols in PcDf.columns:
    CommonDf[Cols.replace('Subs', 'Pc')] = PcDf[Cols]


################# Supplemetal figure 1
CommonDf[CommonDf.HS_Subs == 0].index.to_series().to_csv('HS.null.txt', index=False, header=False)
CommonDf[CommonDf.PT_Subs == 0].index.to_series().to_csv('PT.null.txt', index=False, header=False)
CommonDf[CommonDf.NL_Subs == 0].index.to_series().to_csv('NL.null.txt', index=False, header=False)
CommonDf[CommonDf.GG_Subs == 0].index.to_series().to_csv('GG.null.txt', index=False, header=False)
CommonDf[CommonDf['#1#_Subs'] == 0].index.to_series().to_csv('P1.null.txt', index=False, header=False)


################# 

PvalDf = pd.read_csv('../Subs_gblocks.aln/Simulation/RankTestPvalue.tsv', sep='\t', 
                       index_col = 'Id')

CandPvalDf = PvalDf[(PvalDf[['HS', 'PT', '#1#', 'GG', 'NL']] < 0.01).any(axis=1)]

ExpDf = pd.DataFrame()
ExpDf['HS_Subs'] = PvalDf.Events * 0.1174
ExpDf['PT_Subs'] = PvalDf.Events * 0.1115
ExpDf['#1#_Subs'] = PvalDf.Events * 0.0372
ExpDf['GG_Subs'] = PvalDf.Events * 0.1542
ExpDf['NL_Subs'] = PvalDf.Events * 0.5797

HighExpDf = pd.DataFrame()
n=0
Cols = ExpDf.columns
for Indx, row in CandPvalDf.iterrows():
    HighExpDf.at[Indx, 'Branches'] = ''
    for Col in Cols:
        if CommonDf.at[Indx, Col] <= ExpDf.at[Indx, Col]: continue
        if CandPvalDf.at[Indx, Col.replace('_Subs', '')] >= 0.01: continue
        HighExpDf.at[Indx, 'Branches'] = f"{Col},{HighExpDf.at[Indx, 'Branches']}"
        n += 1
HighExpDf = HighExpDf[HighExpDf.Branches != '']
HighExpDf.Branches.value_counts()

LowExpDf = pd.DataFrame()
n=0
Cols = ExpDf.columns
for Indx, row in CandPvalDf.iterrows():
    LowExpDf.at[Indx, 'Branches'] = ''
    for Col in Cols:
        if CommonDf.at[Indx, Col] >= ExpDf.at[Indx, Col]: continue
        if CandPvalDf.at[Indx, Col.replace('_Subs', '')] >= 0.01: continue
        LowExpDf.at[Indx, 'Branches'] = f"{Col},{LowExpDf.at[Indx, 'Branches']}"
        n += 1
LowExpDf = LowExpDf[LowExpDf.Branches != '']


################# Supplemetal figure 2 Events
with plt.rc_context(dict(sns.axes_style("white", rc = {
    'axes.spines.right': False,
    'axes.spines.top': False,}),
                         **sns.plotting_context("paper", rc={
                             'axes.labelsize':14, 
                             'axes.titlesize': 18,
                             'xtick.labelsize': 12,
                             'ytick.labelsize': 12,}))):    
    plt.close('all')
    plt.clf()
    plt.rcParams['pdf.fonttype'] = 'truetype'

    Data=PvalDf.Events.value_counts().sort_index()
    sns.lineplot(x=Data.index[:100], y=Data.values[:100]).set(
        title='', xlabel = 'Events \'N\'', ylabel = 'Ortholog Families')



CommonDf['lower'] = LowExpDf.Branches.replace(
    ['HS_Subs', 'PT_Subs', '#1#_Subs', 'GG_Subs', 'NL_Subs'], ['Human', 'Chimpanzee','#1#', 'Gorilla','Gibbon'],
    regex=True).replace(',$', '', regex=True)
CommonDf['higher'] = HighExpDf.Branches.replace(
    ['HS_Subs', 'PT_Subs', '#1#_Subs', 'GG_Subs', 'NL_Subs'], ['Human', 'Chimpanzee','#1#', 'Gorilla','Gibbon'],
    regex=True).replace(',$', '', regex=True)


CommonDf.loc[HighExpDf.index.intersection(LowExpDf.index), ['lower', 'higher']]





    
    
#Remove null substituion rows
CommonDf = CommonDf[(CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']] > 0).any(axis=1)]


#Manual curation
CommonDf.loc['ENST00000450565','Comment'] = 'Most divergent gene on the human branch'
CommonDf.loc['ENST00000301633','Comment'] = 'Gibbon misssing starting methionine'
CommonDf.loc['ENST00000391916','Comment'] = 'All substitutions were within a short block'
CommonDf.loc['ENST00000637878', 'Comment'] = 'Among the top 5 divergent genes on the human branch'
CommonDf.loc['ENST00000008938', 'Comment'] = 'Among the top 5 divergent genes on the human branch'
CommonDf.loc['ENST00000259845','Comment'] = 'Among the top 5 divergent genes on the human branch'
CommonDf.loc['ENST00000454136','Comment'] = 'Among the top 5 divergent genes on the human branch'

CommonDf.loc['ENST00000640237','Comment'] = 'First half of the alighment is unrelaible even after filtering'
CommonDf.loc['ENST00000228468','Comment'] = 'First half of the alighment is unrelaible even after filtering'
CommonDf.loc['ENST00000390654','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000265729','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000614943','Comment'] = 'Frame-shift'
CommonDf.loc['ENST00000375464','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000446510','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000373383','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000372398','Comment'] = 'Many substitutions are within a short block'


CommonDf.loc['ENST00000542996','Comment'] = 'Most divergent gene on the chimpanzee branch'
CommonDf.loc['ENST00000296280', 'Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
CommonDf.loc['ENST00000380041', 'Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
CommonDf.loc['ENST00000342995','Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
CommonDf.loc['ENST00000343470','Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'

CommonDf.loc['ENST00000359741','Comment'] = 'Most divergent gene on the #1# branch'
CommonDf.loc['ENST00000397893','Comment'] = 'Many subsituions are within a short block'
CommonDf.loc['ENST00000274520','Comment'] = 'Among the top 5 divergent genes on the #1# branch'
CommonDf.loc['ENST00000598398','Comment'] = 'Short Alignment'
CommonDf.loc['ENST00000370081','Comment'] = 'Short Alignment'
CommonDf.loc['ENST00000299191', 'Comment'] = 'Among the top 5 divergent genes on the #1# branch'
CommonDf.loc['ENST00000292894', 'Comment'] = 'Among the top 5 divergent genes on the #1# branch'
CommonDf.loc['ENST00000456397','Comment'] = 'Many subsituions are within a short block'
CommonDf.loc['ENST00000625099','Comment'] = 'Among the top 5 divergent genes on the #1# branch'


CommonDf.loc['ENST00000370177','Comment'] = 'Short Alignment'
CommonDf.loc['ENST00000255499','Comment'] = 'Most divergent gene on the gorilla branch'
CommonDf.loc['ENST00000359878','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000342790','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000254976', 'Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
CommonDf.loc['ENST00000651546', 'Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
CommonDf.loc['ENST00000649185','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000613760', 'Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
CommonDf.loc['ENST00000263317','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000255977','Comment'] = 'Among the top 5 divergent genes on the gorilla branch'


CommonDf.loc['ENST00000641401','Comment'] = 'Short Alignment'
CommonDf.loc['ENST00000513010','Comment'] = 'Most divergent gene on the gibbon branch'
CommonDf.loc['ENST00000345088', 'Comment'] = 'Among the top 5 divergent genes on the gibbon branch'
CommonDf.loc['ENST00000393330','Comment'] = 'Among the top 5 divergent genes on the gibbon branch'
CommonDf.loc['ENST00000397301','Comment'] = 'Among the top 5 divergent genes on the gibbon branch'
CommonDf.loc['ENST00000523047', 'Comment'] = 'Among the top 5 divergent genes on the gibbon branch'


#save the norm data frame
CommonDf.to_csv('./NormDf.tsv', sep='\t', header=True, index=False)









