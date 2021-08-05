#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Has to be run from the folder where the input files are
"""

import os, pandas as pd, numpy as np
import argparse
from statsmodels.stats import weightstats
from statsmodels.stats.multitest import fdrcorrection


parser = argparse.ArgumentParser(description='Parse the .tsv alignment file to identify substitutions', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f', '--input_csv', required = True, help = 'Full path to the input csv that contians the poisson corrected branch lenghts')
args = vars(parser.parse_args()) #parse arguments 

def FdrCorrection(TheDf):
    """Columnwise fdr correction
    """
    for Col in TheDf.columns:
        TheDf[Col] = fdrcorrection(TheDf[Col])[1]




AvgPcDf = pd.read_csv(args['input_csv'], sep='\t', index_col=0)
BsDf = pd.read_csv('BootStrapSE.tsv', sep='\t', index_col=0)

#change the current working directory
Dr = args['input_csv'].split('/')[-1].split('.')[0]
if os.path.isdir(Dr): pass
else: os.mkdir(Dr)
os.chdir(Dr)

AvgLen = AvgPcDf.mean()
BranchScoreDf = pd.DataFrame(dtype=np.float64)
IdLs = [x.replace('_Subs','') for x in AvgPcDf.columns]

for Indx, row in AvgPcDf.iterrows():
    for Id in IdLs:
        BranchScoreDf.loc[Indx, Id] = abs(row[f"{Id}_Subs"] - AvgLen[f"{Id}_Subs"])#**2
        
BranchScoreDf['TreeFull'] = BranchScoreDf[IdLs].sum(axis=1)
BranchScoreDf.to_csv('BranchScoreDf.tsv', sep='\t', header=True, index=True, index_label='Id')

BraScStatSigDf = pd.DataFrame(dtype=np.float64)
ColLs = list(BranchScoreDf.columns)
for Indx, row in BranchScoreDf.iterrows():
    for Col in ColLs:
        BraScStatSigDf.at[Indx, Col] = weightstats._zstat_generic2(row[Col], BsDf.at[Indx, f"SE_{Col}"], alternative='two-sided',)[-1]/2
BraScStatSigDf.to_csv('BraScStatSigDf.tsv', sep='\t', header=True, index=True, index_label='Id')
FdrCorrection(BraScStatSigDf)
BraScStatSigDf.to_csv('FdrBraScStatSigDf.tsv', sep='\t', header=True, index=True, index_label='Id')