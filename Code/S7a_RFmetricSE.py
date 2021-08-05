#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os, sys, pandas as pd, numpy as np
from statsmodels.stats import weightstats


CommonDf = pd.read_csv('../../CommonDf.tsv', sep='\t', 
                       usecols=['HS', 'overlap', 'HS_Subs', 'PT_Subs', 
                                'GG_Subs', 'NL_Subs', '#1#_Subs',], 
                       index_col = 'HS')
SubCols = [x for x in CommonDf.columns if x.endswith('_Subs')]

#Remove zero alignment families
CommonDf = CommonDf[CommonDf.overlap > 0]
#Remove zero branch-specific substitutions families
CommonDf = CommonDf[(CommonDf[SubCols] > 0).any(axis=1)]

#Poisson correction
PcDf = CommonDf[SubCols].copy(deep=True)
PcDf = PcDf.apply(lambda row: -np.log(1-row/(CommonDf.at[row.name, 'overlap'])), axis=1)
PcDf.to_csv('PoissonCorrectedDf.tsv', sep='\t', header=True, index=True, index_label='Id') 


#Substitution position
SubsPosDf = pd.read_csv('../SubPosDf.tsv', sep='\t', index_col=0, )
SubsPosDf = SubsPosDf.loc[SubsPosDf.index.intersection(PcDf.index),:]
#Fill no substituion branches
SubsPosDf.fillna('-', inplace=True)

def SubsSer(SubsString):
    if SubsString == '-': return([])#empty list
    else:
        SubSer = [int(x) for x in SubsString.split(',')]
        return(SubSer)

StatSignDf = pd.DataFrame(dtype=np.float64)
BsDf = pd.DataFrame(dtype=np.float64)
NormBsDf = pd.DataFrame(dtype=np.float64)
NormHominidBsDf = pd.DataFrame(dtype=np.float64)

m=0
for Indx, row in SubsPosDf.iterrows():
    n=0
    SampleDf = pd.DataFrame()
    Length = int(row['Length'])
    PosSer =  pd.Series(range(1, Length+1, 1))
    HSsubSer = pd.Series(SubsSer(row[ 'HS']))
    PTsubSer = pd.Series(SubsSer(row[ 'PT']))
    P1subSer = pd.Series(SubsSer(row[ '#1#']))
    GGsubSer = pd.Series(SubsSer(row[ 'GG']))
    NLsubSer = pd.Series(SubsSer(row[ 'NL']))
    while n<1000:
        
        SamplePosSer = PosSer.sample(n=Length, replace=True)
        PosValCount = SamplePosSer.value_counts()
        SampleDf.at[n, 'HS'] = PosValCount.loc[PosValCount.index.intersection(HSsubSer.values)].sum()
        SampleDf.at[n, 'PT'] = PosValCount.loc[PosValCount.index.intersection(PTsubSer.values)].sum()
        SampleDf.at[n, '#1#'] = PosValCount.loc[PosValCount.index.intersection(P1subSer.values)].sum()
        SampleDf.at[n, 'GG'] = PosValCount.loc[PosValCount.index.intersection(GGsubSer.values)].sum()
        SampleDf.at[n, 'NL'] = PosValCount.loc[PosValCount.index.intersection(NLsubSer.values)].sum()
        
        n += 1
    Branches = SampleDf.columns
    #Poisson correction
    SampleDf = -np.log(1-SampleDf/Length)
    NormDf = SampleDf.copy(deep=True)
    TreeSum = NormDf.sum(axis=1)   
    SampleDf["TreeLenFull"] = TreeSum
    SampleDf["Var_TreeLenFull"] = (SampleDf["TreeLenFull"] - SampleDf["TreeLenFull"].mean())**2
    BsDf.loc[Indx, "SE_TreeFull"] = ((SampleDf["Var_TreeLenFull"].sum())/999)**0.5
    m += 1

BsDf.to_csv('BootStrapSE.tsv', sep='\t', header=True, index=True, index_label='Id')
    