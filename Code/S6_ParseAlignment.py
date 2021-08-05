#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

from Bio.Alphabet.IUPAC import IUPACProtein
import os, sys, pandas as pd, numpy as np
import argparse

######################  get the input alignment file path
parser = argparse.ArgumentParser(description='Parse the .tsv alignment file to identify substitutions', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f', '--file_path', required = True, help = 'Full path to the aligned file')
parser.add_argument('-g', '--gap', default = False, help = 'If alignment contains gap')
args = vars(parser.parse_args()) #parse arguments 



"""Functions"""
def DrCheck(Dr, make=False): 
    """Check if Dir path exists, if not either exit with error massage or, 
    only if explicittly asked, make the directory"""
    if os.path.isdir(Dr) == False:#not in path
        if make == False:#only check
            sys.exit("\n Directory '{}' is not in the path.".format(Dr))
        else:#make directory
            ParentDr = "/".join(Dr.split("/")[:-1])
            DrCheck(ParentDr)#check if parent exists, bu do not make the parent directory
            sys.stdout.write("""\n Making directory: '{}'""".format(Dr))
            os.mkdir(Dr)

def FlCheck(Dr, Prefix, Suffix, ContinueIfFileExists=False, SendFileLs = True):
    """Check if file exists, return a string"""
    DrCheck(Dr, make=False)#check if directory exists
    FlLs = os.popen(r"""ls {} | grep '{}.*{}$'""".format(Dr, Prefix, Suffix)).read().strip('\n').split('\n')#grep file that ends with the suffix
    FlLs = list(filter(None, FlLs))#check if file exists, filter allows making list of length zero
    if len(FlLs) == 1:# file in directory 
        if ContinueIfFileExists == True: return('continue')
        else: return(Dr+"/"+FlLs[0])#return file name
    elif len(FlLs) == 0:# file not in directory 
        sys.exit("""\n File not found.\n
                     Directory : {} \n
                     Prefix: {} \n
                     Suffix: {} \n""".format(Dr, Prefix, Suffix))
    elif len(FlLs) > 1:#multiple files, exit
        if SendFileLs == True: return(FlLs)
        sys.exit("\n File identifier not unique.\nFiles; {}".format(FlLs))

def CreateSubsString(SubsDf, SubsColLs):
    """Given a dataframe of substituted sites and substituted columns list.
    Returns substitution string
    """
    if len(SubsDf) == 0: return(np.nan)
    Df = SubsDf.copy(deep=True) 
    ConsCols = [x for x in Df.columns if x not in SubsColLs]
    Df['Pos'] =  Df.index + 1
    Df['Subs'] = Df.apply(lambda row: f"{row['Pos']}:{row[SubsColLs[0]]}|{row[ConsCols[0]]}",axis=1)
    return((",").join(Df.Subs.values))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

AlnFl = args['file_path']
InDir = "/".join(AlnFl.split("/")[:-1])
AlnFlName = AlnFl.split("/")[-1]
Prefix = ".".join(AlnFlName.split('.')[:-1])
Suffix = AlnFlName.split('.')[-1]
FlCheck(InDir, Prefix, Suffix)
OutDir = f"{InDir}/Subs_{Prefix}"
DrCheck(OutDir, make = True)


GblockDf = pd.read_csv(AlnFl, sep='\t')
AlignedCols = [x for x in GblockDf.columns if x.endswith('_al')]
IdCols = [x for x in GblockDf.columns if not x.endswith('_al')]
SubResDf = pd.DataFrame(dtype=np.str, columns=["Length"] + IdCols)
IdenityDf = pd.DataFrame(dtype=np.int, columns=list(IUPACProtein.letters)) #list(IUPACProtein.letters) + ['X', '*']
TransIdSp = GblockDf.columns[0]
for Indx,row in GblockDf.iterrows():
    TransId = row[0].replace('>','')
    print(TransId)
    Alignment = row[AlignedCols]
    AlignDf = pd.DataFrame(np.array([list(rec) for rec in Alignment]), index=Alignment.index).transpose()
    if args['gap'] != False:
        AlignDf = AlignDf.dropna(axis=1, how='any')
        AlignDf = AlignDf[~(AlignDf == '-').any(axis=1)]
        AlignDf = AlignDf[~(AlignDf == 'X').any(axis=1)]
        AlignDf.reset_index(drop=True, inplace=True)
    UnqRowVal = AlignDf.nunique(axis=1)
    TwoResDf = AlignDf.loc[UnqRowVal[UnqRowVal == 2].index, :].copy(deep=True)
    PhyloConsDf = TwoResDf[~((TwoResDf.HS_al ==  TwoResDf.GG_al) & (TwoResDf.PT_al ==  TwoResDf.NL_al))]
    PhyloConsDf = PhyloConsDf[~((PhyloConsDf.HS_al ==  PhyloConsDf.NL_al) & (PhyloConsDf.PT_al ==  PhyloConsDf.GG_al))]
    SubsP1 = PhyloConsDf[((PhyloConsDf.HS_al ==  PhyloConsDf.PT_al) & (PhyloConsDf.NL_al ==  PhyloConsDf.GG_al))]
    SubsHS = PhyloConsDf[((PhyloConsDf.HS_al !=  PhyloConsDf.PT_al) & (PhyloConsDf.PT_al ==  PhyloConsDf.GG_al))]
    SubsPT = PhyloConsDf[((PhyloConsDf.HS_al !=  PhyloConsDf.PT_al) & (PhyloConsDf.HS_al ==  PhyloConsDf.GG_al))]
    SubsGG = PhyloConsDf[((PhyloConsDf.HS_al ==  PhyloConsDf.NL_al) & (PhyloConsDf.NL_al !=  PhyloConsDf.GG_al))]
    SubsNL = PhyloConsDf[((PhyloConsDf.HS_al ==  PhyloConsDf.GG_al) & (PhyloConsDf.NL_al !=  PhyloConsDf.GG_al))]
    
    SubResDf.loc[TransId, 'Length'] = np.str(len(AlignDf.index))
    SubResDf.loc[TransId, 'HS'] = CreateSubsString(SubsHS, ['HS_al',])
    SubResDf.loc[TransId, 'PT'] = CreateSubsString(SubsPT, ['PT_al',])
    SubResDf.loc[TransId, '#1#'] = CreateSubsString(SubsP1, ['HS_al', 'PT_al'])
    SubResDf.loc[TransId, 'GG'] = CreateSubsString(SubsGG, ['GG_al',])
    SubResDf.loc[TransId, 'NL'] = CreateSubsString(SubsNL, ['NL_al',])
    

SubResDf.to_csv("{}/SubResDf.tsv".format(OutDir), sep='\t', header=True, index=True, index_label='Id')
os.popen(f"""sed 's/:[A-Z]|[A-Z]//g' {OutDir}/SubResDf.tsv > {OutDir}/SubPosDf.tsv
         rm {OutDir}/SubResDf.tsv
         """).read()

