#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 15:31:47 2021

@author: prabh
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
import numpy as np
import os
from statsmodels.stats.multitest import multipletests



#os.chdir('/groups/prodiv/projects/2019-10-21_mammals/data/msa/Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys/Homo_sapiens/Subs_gblocks.aln/Simulation')

PvalDf = pd.read_csv('RankTestPvalue.tsv', sep='\t', 
                       index_col = 'Id')

CommonDf = pd.read_csv('../../CommonDf.tsv', sep='\t', index_col=0)
CommonDf.dropna(inplace=True, subset=['overlap']) #12618
AaLs = [x for x in list(CommonDf.columns) if x.endswith('_aa')]


CommonDf['Min'] = CommonDf[AaLs].apply(lambda row: min(row), axis=1)
CommonDf['Max'] = CommonDf[AaLs].apply(lambda row: max(row), axis=1)
################# figure 2
from scipy.stats import linregress
linregress(CommonDf.Min, CommonDf.Max)
#LinregressResult(slope=1.0263976539355473, intercept=23.951985096258454, rvalue=0.9741702021024219, pvalue=0.0, stderr=0.0021182283172350605)

CommonDf['MaxDif'] = CommonDf[AaLs].apply(lambda row: (max(row) - min(row)), axis=1)
CommonDf['RelMaxDif'] = CommonDf['MaxDif']*100/CommonDf['Min']

CommonDf['RelOverlap'] = CommonDf['overlap']*100/CommonDf['Min']
CommonDf['%MinId'] = CommonDf['AbsId']*100/CommonDf['Min']


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
LowExpDf.Branches.value_counts()
"""
HighExpDf.Branches.value_counts()
Out[60]: 
GG_Subs,             137
NL_Subs,             121
PT_Subs,             101
HS_Subs,              92
#1#_Subs,             79
PT_Subs,HS_Subs,       1
GG_Subs,#1#_Subs,      1
Name: Branches, dtype: int64

LowExpDf.Branches.value_counts()
Out[61]: 
NL_Subs,                     117
HS_Subs,                      17
GG_Subs,                      15
PT_Subs,                       8
#1#_Subs,                      2
NL_Subs,GG_Subs,               2
GG_Subs,HS_Subs,               1
GG_Subs,#1#_Subs,HS_Subs,      1
Name: Branches, dtype: int64
"""

#Both High and Low
len(HighExpDf.index.intersection(LowExpDf.index))
#One sp overlap
LowInOneBr = LowExpDf[[x in ['HS_Subs,', 'PT_Subs,', '#1#_Subs,', 'GG_Subs,', 'NL_Subs,'] for x in LowExpDf.Branches]]
HighInOneBr = HighExpDf[[x in ['HS_Subs,', 'PT_Subs,', '#1#_Subs,', 'GG_Subs,', 'NL_Subs,'] for x in HighExpDf.Branches]]
HighInOneBr.index.union(LowInOneBr.index)
HighInOneBr.index.intersection(LowInOneBr.index)
LowHighDf =  LowInOneBr.loc[HighInOneBr.index.intersection(LowInOneBr.index),:].copy()
LowHighDf = LowHighDf + HighInOneBr.loc[HighInOneBr.index.intersection(LowInOneBr.index),:]
LowHighDf.Branches.value_counts()

#plots
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

Data=PvalDf.loc[HighInOneBr.index,'Events'].value_counts().sort_index()
sns.lineplot(x=Data.index[Data.index <= 100], y=Data.values[Data.index <= 100]).set(
    title='High Than Expected', xlabel = 'Events \'N\'', ylabel = 'Ortholog Families')

Data=PvalDf.loc[LowInOneBr.index,'Events'].value_counts().sort_index()
sns.lineplot(x=Data.index[Data.index <= 100], y=Data.values[Data.index <= 100]).set(
    title='Low Than Expected', xlabel = 'Events \'N\'', ylabel = 'Ortholog Families')


CommonDf['lower'] = LowExpDf.Branches.replace(
    ['HS_Subs', 'PT_Subs', '#1#_Subs', 'GG_Subs', 'NL_Subs'], ['Human', 'Chimpanzee','#1#', 'Gorilla','Gibbon'],
    regex=True).replace(',$', '', regex=True)
CommonDf['higher'] = HighExpDf.Branches.replace(
    ['HS_Subs', 'PT_Subs', '#1#_Subs', 'GG_Subs', 'NL_Subs'], ['Human', 'Chimpanzee','#1#', 'Gorilla','Gibbon'],
    regex=True).replace(',$', '', regex=True)


CommonDf.loc[HighExpDf.index.intersection(LowExpDf.index), ['lower', 'higher']]
CommonDf.loc[HighExpDf.index.intersection(LowExpDf.index), ['lower', 'higher']].value_counts()


if os.path.isdir('candidates') == False: os.mkdir('candidates')
CommonDf.loc[HighExpDf[HighExpDf.Branches == 'HS_Subs,'].index, ['Gene', 'Description', '%HS_Subs', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%HS_Subs',ascending=False).to_csv('candidates/HighExp_Human.tsv', sep='\t', header=['Gene', 'Description', 'Branch-specific % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[HighExpDf[HighExpDf.Branches == 'PT_Subs,'].index, ['PT', 'Gene', 'Description', '%PT_Subs', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%PT_Subs',ascending=False).to_csv('candidates/HighExp_Chimpanzee.tsv', sep='\t', header=['Chimpanzee transcript', 'Gene', 'Description', 'Branch-specific % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[HighExpDf[HighExpDf.Branches == '#1#_Subs,'].index, ['Gene', 'Description', '%#1#_Subs',  'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%#1#_Subs',ascending=False).to_csv('candidates/HighExp_P1.tsv', sep='\t', header=['Gene', 'Description', 'Branch-specific % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[HighExpDf[HighExpDf.Branches == 'GG_Subs,'].index, ['GG', 'Gene', 'Description', '%GG_Subs', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%GG_Subs',ascending=False).to_csv('candidates/HighExp_Gorilla.tsv', sep='\t', header=['Gorilla transcript', 'Gene', 'Description', 'Branch-specific % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[HighExpDf[HighExpDf.Branches == 'NL_Subs,'].index, ['NL', 'Gene', 'Description', '%NL_Subs', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%NL_Subs',ascending=False).to_csv('candidates/HighExp_Gibbon.tsv', sep='\t', header=['Gibbon transcript', 'Gene', 'Description', 'Branch-specific % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')

MixedHigh = HighExpDf[[x not in ['HS_Subs,', 'PT_Subs,', 'GG_Subs,', '#1#_Subs,', 'NL_Subs,'] for x in HighExpDf.Branches]]
CommonDf.loc[MixedHigh.index, ['Gene', 'Description', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower', 'higher']].to_csv('candidates/MixedHighExp.tsv', sep='\t', header=['Gene', 'Description','Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')

CommonDf.loc[LowExpDf[LowExpDf.Branches == 'HS_Subs,'].index, ['Gene', 'Description', '%HS_Subs', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= '%HS_Subs',ascending=True).to_csv('candidates/LowExp_Human.tsv', sep='\t', header=['Gene', 'Description', 'Human % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[LowExpDf[LowExpDf.Branches == 'PT_Subs,'].index, ['PT', 'Gene', 'Description', '%PT_Subs',  'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= '%PT_Subs',ascending=True).to_csv('candidates/LowExp_Chimpanzee.tsv', sep='\t', header=['Chimpanzee transcript', 'Gene', 'Description', 'Chimp % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[LowExpDf[LowExpDf.Branches == '#1#_Subs,'].index, ['PT', 'Gene', 'Description', '%#1#_Subs',  'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= '%#1#_Subs',ascending=True).to_csv('candidates/LowExp_P1.tsv', sep='\t', header=['Chimpanzee transcript', 'Gene', 'Description', 'Chimp % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[LowExpDf[LowExpDf.Branches == 'GG_Subs,'].index, ['GG', 'Gene', 'Description', '%GG_Subs',  'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= '%GG_Subs',ascending=True).to_csv('candidates/LowExp_Gorilla.tsv', sep='\t', header=['Gorilla transcript', 'Gene', 'Description', 'Gorilla % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[LowExpDf[LowExpDf.Branches == 'NL_Subs,'].index, ['NL', 'Gene', 'Description', '%NL_Subs',  'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= '%NL_Subs',ascending=True).to_csv('candidates/LowExp_Gibbon.tsv', sep='\t', header=['Gibbon transcript', 'Gene', 'Description', 'Gibbon % subs per site', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')


#write HTML alignment files
if os.path.isdir('HTML') == False: os.mkdir('HTML')
for Id in CandPvalDf.index:
    GeneName = CommonDf.at[Id, 'Gene']
    os.popen(f"""python3 /groups/prodiv/projects/2019-10-21_mammals/code/2020.12.01/Supli_S5_ExtractSingleAlignmentFromGblocksHTML.py -i {Id} -f ../../Gblocks.html -o HTML/{GeneName}___""").read()

os.popen(r"rm HTML/*.bak").read()





