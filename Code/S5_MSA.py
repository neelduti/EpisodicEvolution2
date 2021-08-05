#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input:
    
    1. Get the speceies names and the focal species name.
    2. Go to their synteny folder.
    3. Check the id.<Focal Species>.txt file for identifiers, ensure that the first line contains the focal species.
    4. Get the suffix of pairwise files with the focal species(collinearity or best blast match).
    5. Read these files as DataFrames and drop focal species duplicates.
    6. Make a combined dataframe only with focal genes with match in every comparision.
    7. 

"""
import time
import os, sys, pandas as pd, numpy as np
from Bio import AlignIO
import argparse
from Bio import SeqIO


######################  get the input species names and directory #######################


parser = argparse.ArgumentParser(description='Parses pairwise species orthologs, finds orthologs across species and does MSA', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-s', '--species', required = True, help = 'Name of the species, without space and seperated by comma, in the same order as it exists in the synteny directory')
parser.add_argument('-S', '--syn_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/synteny', help = 'folder that has blastp files')
parser.add_argument('-f', '--focal', required = True, help = 'Name of focal sp')
parser.add_argument('-p', '--pairwise_suf', default = 'CollinearDf.tsv')
parser.add_argument('-I', '--sp_id_file', default = 'Id.txt')
parser.add_argument('-F', '--fasta_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/primTran', help = 'folder containing fasta files')
parser.add_argument('-P', '--pep_suf', default = 'peptide.fa')
parser.add_argument('-o', '--out_fol', default = '/groups/prodiv/projects/2019-10-21_mammals/data/msa', help = 'folder that will have the output')
parser.add_argument('-a', '--ancestral', default = 0)
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

def DropDuplicates(FlName, FocalId, OthSpId):
    """Sort the gene pairs by the e-value (pairwise blast hit) and the length of their syntenic 
    blocks and keep only the first instance of each focal species gene.
    """
    PwDf = pd.read_csv(FlName, sep='\t')
    #sort values according to e_value of the alignment which is based on the number of genes in the alignment
    PwDf = PwDf.sort_values(['e_value', 'IdRes', 'rank'], ascending=[True, False, True], axis=0).copy()
    #keep the first match which occurs in the largest syntenic block
    PwDf.drop_duplicates(subset = FocalId, keep = 'first', inplace=True)
    PwDf.drop_duplicates(subset = OthSpId, keep = 'first', inplace=True)
    PwDf.set_index(FocalId, inplace=True)
    return(PwDf)

def WriteToAllAlignmentDf(indx, alignDf, SpIdLs, FastaIdLs, AllAlignmntDf):
    """Given the Alignment dataframe, species id and output dataframe,
    it writes the alignment to the output dataframe.    
    """
    n=0
    for Id in FastaIdLs:
        AllAlignmntDf.loc[indx, SpIdLs[n]] = '>' + Id
        AllAlignmntDf.loc[indx, '{}_al'.format(SpIdLs[n])] = ''.join(alignDf[n])
        n += 1

def MSA(CommonDf):
    """Create MSA for each row of the CommonDf, return alignment identity as 
    well as the alignment itself as tsv files.
    """
    Sp1Id = CommonDf.index.name
    Sp2Id = SpIdDf.id[1]
    Sp3Id = SpIdDf.id[2]
    Sp4Id = SpIdDf.id[3]
    SpIdLs = [Sp1Id, Sp2Id, Sp3Id, Sp4Id]
    AncIdLs = [Sp1Id,'#1#',Sp2Id, '#2#',Sp3Id, '#3#',Sp4Id]
    PepFlSp1 = FlCheck(Dr = FastaDr, Prefix = SpNameDf.at[Sp1Id, 'species'], Suffix = args["pep_suf"])
    PepFlSp2 = FlCheck(Dr = FastaDr, Prefix = SpNameDf.at[Sp2Id, 'species'], Suffix = args["pep_suf"])
    PepFlSp3 = FlCheck(Dr = FastaDr, Prefix = SpNameDf.at[Sp3Id, 'species'], Suffix = args["pep_suf"])
    PepFlSp4 = FlCheck(Dr = FastaDr, Prefix = SpNameDf.at[Sp4Id, 'species'], Suffix = args["pep_suf"])
    PepSp1 = SeqIO.index(PepFlSp1, 'fasta')
    PepSp2 = SeqIO.index(PepFlSp2, 'fasta')
    PepSp3 = SeqIO.index(PepFlSp3, 'fasta')
    PepSp4 = SeqIO.index(PepFlSp4, 'fasta')
    aaSp1 = '{}_aa'.format(Sp1Id)
    aaSp2 = '{}_aa'.format(Sp2Id)
    aaSp3 = '{}_aa'.format(Sp3Id)
    aaSp4 = '{}_aa'.format(Sp4Id)

    AllMafftAlnDf = pd.DataFrame()
    AllGblocksAlnDf = pd.DataFrame()
    # AllPrankAncAlnDf = pd.DataFrame()
    # AncSubResDf = pd.DataFrame(columns=AncIdLs)
    # AncSubCountDf = pd.DataFrame(columns=AncIdLs)
    SubResDf = pd.DataFrame()


    
    #remove if the combined html file already exists
    ConcatHTMLfl = "{}/Gblocks.html".format(OutFocalSpDataDir)
    os.popen(r"echo -n > {}".format(ConcatHTMLfl)).read()
    
    #remove if combined ancestral events file exists
    # ConcatAncEventsFl = "{}/AncEvents.txt".format(OutFocalSpDataDir)
    # os.popen(r"echo -n > {}".format(ConcatAncEventsFl)).read()
    
    
    n = 0
    for indx, row in CommonDf.iterrows():
        os.popen(r"rm -rf tmp*").read()
        with open('tmp.tree.nwk', 'w') as fh:
            fh.write(r"((({}, {}), {}), {});".format(indx, row[Sp2Id], row[Sp3Id], row[Sp4Id]))
        with open('tmp.protein.fa', 'w') as fh: 
            #make fasta file for alignment
            SeqIO.write(PepSp1[indx], fh, 'fasta')
            SeqIO.write(PepSp2[row[Sp2Id]], fh, 'fasta')
            SeqIO.write(PepSp3[row[Sp3Id]], fh, 'fasta')
            SeqIO.write(PepSp4[row[Sp4Id]], fh, 'fasta')
            
        #write sequence length in CommonDf
        if PepSp1[indx].seq[-1] == '*': CommonDf.at[indx,aaSp1] = len(PepSp1[indx].seq)-1
        else: CommonDf.at[indx,aaSp1] = len(PepSp1[indx].seq)
        if PepSp2[row[Sp2Id]].seq[-1] == '*': CommonDf.at[indx,aaSp2] = len(PepSp2[row[Sp2Id]].seq)-1
        else: CommonDf.at[indx,aaSp2] = len(PepSp2[row[Sp2Id]].seq)
        if PepSp3[row[Sp3Id]].seq[-1] == '*': CommonDf.at[indx,aaSp3] = len(PepSp3[row[Sp3Id]].seq)-1
        else: CommonDf.at[indx,aaSp3] = len(PepSp3[row[Sp3Id]].seq)
        if PepSp4[row[Sp4Id]].seq[-1] == '*': CommonDf.at[indx,aaSp4] = len(PepSp4[row[Sp4Id]].seq)-1
        else: CommonDf.at[indx,aaSp4] = len(PepSp4[row[Sp4Id]].seq)
        
        #run mafft on the peptide file, fix gap open penalty
        os.popen(r"""    
        mafft --quiet --op 3 tmp.protein.fa > tmp.aln.mafft.fa
        sed -i 's/ <unknown description>//g' tmp.aln.mafft.fa
        """).read()

        """alignment parsing 1"""
        alignmnt = AlignIO.read('tmp.aln.mafft.fa', 'fasta')
        alignDf = pd.DataFrame(np.array([list(rec) for rec in alignmnt])).transpose()
        #write to alignment file
        FastaIdLs = [x.id for x in alignmnt]
        WriteToAllAlignmentDf(indx, alignDf,  SpIdLs, FastaIdLs, AllMafftAlnDf)
        
        #run gblocks, remove gap columns and blocks with ....
        os.popen(r"""    
        Gblocks tmp.aln.mafft.fa -b2=4 -b1=4 -b3=2 -b5=n > tmp.gb.fa.out
        """).read()
        
        """alignment parsing 2"""
        #check if gblocks saves atleast one position in the tmp.gb.fa.out file
        GblocksAlnLs = os.popen(r"""grep 'Gblocks alignment:' tmp.gb.fa.out""").read().split()
        if int(GblocksAlnLs[2]) <= 0: continue 
    
        alignmnt = AlignIO.read('tmp.aln.mafft.fa-gb', 'fasta')
        GBalignDf = pd.DataFrame(np.array([list(rec) for rec in alignmnt])).transpose()
        #write to alignment file
        FastaIdLs = [x.id for x in alignmnt]
        WriteToAllAlignmentDf(indx, GBalignDf,  SpIdLs, FastaIdLs, AllGblocksAlnDf)
        
        #Append Gblocks html files into one
        os.popen(r"""
                 sed -i 's/Gblocks 0.91b Results/{}/g' tmp.aln.mafft.fa-gb.htm
                 cat tmp.aln.mafft.fa-gb.htm >> {}
                 """.format(indx, ConcatHTMLfl)).read()
        
        # #run prank to extract ancestral events (-keep)
        # if int(args['ancestral']) == 1: 
        #     os.popen(r"""    
        #              prank -d=tmp.aln.mafft.fa-gb -t=tmp.tree.nwk -o=tmp.prank.aln.fa -showanc -showevents -once -keep
        #              """).read()
        # #check if ancestral events file is generate, this file is not generated if there are no sunstituions to report
        # if int(args['ancestral']) == 1 and os.popen(r"ls tmp.prank.aln.fa.events").read() != '':   
        #     alignmnt = AlignIO.read('tmp.prank.aln.fa.anc.fas', 'fasta')
        #     alignDf = pd.DataFrame(np.array([list(rec) for rec in alignmnt])).transpose()
        #     #write to alignment file
        #     FastaIdLs = [x.id for x in alignmnt]
        #     WriteToAllAlignmentDf(indx, alignDf,  AncIdLs, FastaIdLs, AllPrankAncAlnDf)
        
        #     """Count and store substiotion on each branch and make substitution matrix for each branch
        #     """
            
        #     #make a single text file out of all ancestral events
        #     os.popen(r"""
        #               cat tmp.prank.aln.fa.events >> {}
        #               """.format(ConcatAncEventsFl)).read()
            
        #     #store ancestral alignment in a dataframe
        #     SeqDf = pd.DataFrame({'SpId':AncIdLs, 'FastaId':FastaIdLs}) #create fastaid, species dataframe
        #     SeqDf.set_index('FastaId', inplace=True) #make fasta id the index column
        #     #read the ancestral events file as a pandas series
        #     AncEvnetsSer = pd.Series(os.popen(r"cat tmp.prank.aln.fa.events").read().split('\n'))
        #     AncEvnetsSer = AncEvnetsSer[7:]
        #     AncEvnetsSer = AncEvnetsSer[~(AncEvnetsSer == '')] # drop empty lines
        #     for line in AncEvnetsSer:
        #         line = line.split()
        #         #print(line)
        #         if line[0] == 'branch':
        #             Branch = SeqDf.at[line[1], 'SpId']
        #             # print(Branch)
        #             AncSubResDf.loc[indx, Branch] = ''
        #             AncSubCountDf.at[indx, Branch] = 0
        #         if len(line)== 4 and line[2] == '->':
        #             AncSubResDf.at[indx, Branch] = AncSubResDf.at[indx, Branch] + ("{}:{}|{},".format(line[0], line[1], line[3]))
        #             AncSubCountDf.at[indx, Branch] = AncSubCountDf.at[indx, Branch] + 1

        
        #overlap length
        CommonDf.at[indx, 'overlap'] = len(GBalignDf)
        UnqRowVal = GBalignDf.nunique(axis=1)
        #absolute identity
        AbsIdenDf = GBalignDf.loc[UnqRowVal[UnqRowVal == 1].index, :]#.copy(deep=True)
        RemainDf = GBalignDf.loc[UnqRowVal[UnqRowVal > 1].index, :].copy(deep=True)      
        #Subs
        UnqRowValRemainDf = RemainDf.nunique(axis=1)
        RemainDf.replace(to_replace=np.nan, value='-', inplace=True)
        SubsDf = RemainDf[(UnqRowValRemainDf > 1) & (UnqRowValRemainDf < 4)].copy(deep=True)
        SubResDf.at[indx, 'Subs'] = ''
        for pos, residues in SubsDf.iterrows():
            SubResDf.at[indx, 'Subs'] = SubResDf.at[indx, 'Subs'] + ("{}:{},".format(pos + 1, '|'.join(residues)))
            
        """single substituion event"""
        SingleSubDf = RemainDf[UnqRowValRemainDf == 2].copy(deep=True)
        #Subs clearly at the terminal branch
        SubInSp1Df =  SingleSubDf[(SingleSubDf[0] != SingleSubDf[1]) & (SingleSubDf[1] == SingleSubDf[2]) & (SingleSubDf[2] == SingleSubDf[3])]#.copy(deep=True)
        SubInSp2Df =  SingleSubDf[(SingleSubDf[0] != SingleSubDf[1]) & (SingleSubDf[0] == SingleSubDf[2]) & (SingleSubDf[2] == SingleSubDf[3])]#.copy(deep=True)
        SubInSp3Df =  SingleSubDf[(SingleSubDf[2] != SingleSubDf[1]) & (SingleSubDf[1] == SingleSubDf[0]) & (SingleSubDf[0] == SingleSubDf[3])]#.copy(deep=True)
        SubInSp4Df =  SingleSubDf[(SingleSubDf[3] != SingleSubDf[1]) & (SingleSubDf[1] == SingleSubDf[2]) & (SingleSubDf[2] == SingleSubDf[0])]#.copy(deep=True)

        #derived subs or convergent substitution at the terminal branch
        SubInAnc1Df  = SingleSubDf[(SingleSubDf[0] == SingleSubDf[1]) & (SingleSubDf[2] == SingleSubDf[3])]#.copy(deep=True)
        ConvergentSubDf = SingleSubDf.loc[SingleSubDf.index.difference(SubInSp1Df.index.union(SubInSp2Df.index.union(SubInSp3Df.index.union(SubInSp4Df.index.union(SubInAnc1Df.index))))), :]#.copy(deep=True))
        """Multiple substitutions"""
        #only a single pair has identity
        ThreeRes = RemainDf[UnqRowValRemainDf == 3]#.copy(deep=True)
        OnlyInGpId = ThreeRes[(ThreeRes[0] == ThreeRes[1])]#.copy(deep=True)
        OnlyOutGpId = ThreeRes[(ThreeRes[2] == ThreeRes[3])]#.copy(deep=True)
        OneInOutId = ThreeRes.loc[ThreeRes.index.difference(OnlyInGpId.index.union(OnlyOutGpId.index)), :]#.copy(deep=True)
        #NoId
        NoIdDf = RemainDf[UnqRowValRemainDf == 4]#.copy(deep=True)
        
        CommonDf.at[indx,'AbsId'] = len(AbsIdenDf)
        CommonDf.at[indx,'Subs'] = len(SubsDf)
        CommonDf.at[indx,'{}_Subs'.format(Sp1Id)] = len(SubInSp1Df)
        CommonDf.at[indx,'{}_Subs'.format(Sp2Id)] = len(SubInSp2Df)
        CommonDf.at[indx,'{}_Subs'.format(Sp3Id)] = len(SubInSp3Df)
        CommonDf.at[indx,'{}_Subs'.format(Sp4Id)] = len(SubInSp4Df)
        CommonDf.at[indx,'#1#_Subs'.format(Sp4Id)] = len(SubInAnc1Df)
        CommonDf.at[indx,'Convergent_Subs'.format(Sp4Id)] = len(ConvergentSubDf)
        CommonDf.at[indx,'OnlyInGpId'] = len(OnlyInGpId)
        CommonDf.at[indx,'OnlyOutGpId'] = len(OnlyOutGpId)
        CommonDf.at[indx,'OneInOutId'] = len(OneInOutId)
        CommonDf.at[indx,'NoId'] = len(NoIdDf)

        
        """Pairwise alignment parsing"""
        # CommonDf.at[indx, '{}_{}_Id'.format(Sp1Id, Sp2Id)
        #             ] = PairwiseAlignmentParsing(GBalignDf, 0, 1)
        # CommonDf.at[indx, '{}_{}_Id'.format(Sp1Id, Sp3Id)
        #             ] = PairwiseAlignmentParsing(GBalignDf, 0, 2)
        # CommonDf.at[indx, '{}_{}_Id'.format(Sp1Id, Sp4Id)
        #             ] = PairwiseAlignmentParsing(GBalignDf, 0, 3)
        # CommonDf.at[indx, '{}_{}_Id'.format(Sp2Id, Sp3Id)
        #             ] = PairwiseAlignmentParsing(GBalignDf, 1, 2)
        # CommonDf.at[indx, '{}_{}_Id'.format(Sp2Id, Sp4Id)
        #             ] = PairwiseAlignmentParsing(GBalignDf, 1, 3)
        # CommonDf.at[indx, '{}_{}_Id'.format(Sp3Id, Sp4Id)
        #             ] = PairwiseAlignmentParsing(GBalignDf, 2, 3)

        # n += 1
        # if n > 10: break

    
    CommonDf['%AbsId'] = CommonDf['AbsId']*100/CommonDf['overlap']
    CommonDf['%Subs'] = 100 - CommonDf['%AbsId']
    CommonDf['%{}_Subs'.format(Sp1Id)] = CommonDf['{}_Subs'.format(Sp1Id)]*100/CommonDf['overlap']
    CommonDf['%{}_Subs'.format(Sp2Id)] = CommonDf['{}_Subs'.format(Sp2Id)]*100/CommonDf['overlap']
    CommonDf['%{}_Subs'.format(Sp3Id)] =  CommonDf['{}_Subs'.format(Sp3Id)]*100/CommonDf['overlap']
    CommonDf['%{}_Subs'.format(Sp4Id)] = CommonDf['{}_Subs'.format(Sp4Id)]*100/CommonDf['overlap']
    CommonDf['%#1#_Subs'] = CommonDf['#1#_Subs']*100/CommonDf['overlap']
    CommonDf['%Convergent_Subs'] = CommonDf['Convergent_Subs']*100/CommonDf['overlap']
    CommonDf['%OnlyInGpId'] = CommonDf['OnlyInGpId']*100/CommonDf['overlap']
    CommonDf['%OnlyOutGpId'] = CommonDf['OnlyOutGpId']*100/CommonDf['overlap']
    CommonDf['%OneInOutId'] = CommonDf['OneInOutId']*100/CommonDf['overlap']
    CommonDf['%NoId'] = CommonDf['NoId']*100/CommonDf['overlap']
    # CommonDf['%_{}_{}_Id'.format(Sp1Id, Sp2Id)] = CommonDf['{}_{}_Id'.format(Sp1Id, Sp2Id)]*100/CommonDf['overlap']
    # CommonDf['%_{}_{}_Id'.format(Sp1Id, Sp3Id)] = CommonDf['{}_{}_Id'.format(Sp1Id, Sp3Id)]*100/CommonDf['overlap']
    # CommonDf['%_{}_{}_Id'.format(Sp1Id, Sp4Id)] = CommonDf['{}_{}_Id'.format(Sp1Id, Sp4Id)]*100/CommonDf['overlap']
    # CommonDf['%_{}_{}_Id'.format(Sp2Id, Sp3Id)] = CommonDf['{}_{}_Id'.format(Sp2Id, Sp3Id)]*100/CommonDf['overlap']
    # CommonDf['%_{}_{}_Id'.format(Sp2Id, Sp4Id)] = CommonDf['{}_{}_Id'.format(Sp2Id, Sp4Id)]*100/CommonDf['overlap']
    # CommonDf['%_{}_{}_Id'.format(Sp3Id, Sp4Id)] = CommonDf['{}_{}_Id'.format(Sp3Id, Sp4Id)]*100/CommonDf['overlap']
    CommonDf.to_csv("{}/CommonDf.tsv".format(OutFocalSpDataDir), sep='\t', header=True, index=True)
    AllMafftAlnDf.to_csv("{}/mafft.aln.tsv".format(OutFocalSpDataDir), sep='\t', header=True, index=False)
    AllGblocksAlnDf.to_csv("{}/gblocks.aln.tsv".format(OutFocalSpDataDir), sep='\t', header=True, index=False)
    # if int(args['ancestral']) == 1: AllPrankAncAlnDf.to_csv("{}/prank.anc.aln.tsv".format(OutFocalSpDataDir), sep='\t', header=True, index=False)
    # AncSubResDf.to_csv("{}/AncStateSubsRes.tsv".format(OutFocalSpDataDir), sep='\t', header=True, index=True, na_rep='-')
    # AncSubCountDf.to_csv("{}/AncStateSubsCount.tsv".format(OutFocalSpDataDir), sep='\t', header=True, index=True, na_rep=0)   


def PairwiseAlignmentParsing(alignDf, Sp1Col, Sp2Col):
    """Parse pairwise alignment from the MSA returns identity and substitution"""
    PwDf = alignDf[[Sp1Col, Sp2Col]]
    Id = len(PwDf[PwDf[Sp1Col] == PwDf[Sp2Col]])
    return(Id)
    
    




"""Input"""
SpList = args["species"].split(',')
if len(SpList) < 2: sys.exit("""\n Less than two species given as input """)
FocalSp = args["focal"]
SynDr = args["syn_dir"]
DrCheck(SynDr)
AllSpStr = args["species"].replace(' ', '_')
SpSynDir = SynDr + "/" + AllSpStr
DrCheck(SpSynDir)
FastaDr = args["fasta_dir"]
DrCheck(FastaDr)
OutDataDr = args["out_fol"]
DrCheck(OutDataDr)
OutDocDr = OutDataDr.replace('data', 'doc')
DrCheck(OutDocDr)

"""Output"""
OutMsaDataDir = OutDataDr + "/" + AllSpStr
DrCheck(OutMsaDataDir, make = True)
OutFocalSpDataDir = OutMsaDataDir + "/" + FocalSp
DrCheck(OutFocalSpDataDir, make = True)
OutMsaDocDir = OutDocDr + "/" + AllSpStr
DrCheck(OutMsaDocDir, make = True)
OutFocalSpDocDir = OutMsaDocDir + "/" + FocalSp
DrCheck(OutFocalSpDocDir, make = True)
os.chdir(OutFocalSpDocDir)

"""Read SpIdFl as a dataframe
"""
SpIdFl = FlCheck(Dr=SpSynDir, Prefix = args["sp_id_file"], Suffix = '')
SpIdDf = pd.read_csv(SpIdFl, sep='\t', header = None, index_col = 1, names=('id', 'species'))
SpNameDf = pd.read_csv(SpIdFl, sep='\t', header = None, index_col = 0, names=('id', 'species'))
FocalId = SpIdDf.at[FocalSp, 'id']
#Pairwise file list
PwFlLs = FlCheck(Dr = SpSynDir, Prefix = FocalId, Suffix = args['pairwise_suf'], SendFileLs = True)
# sys.stdout.write(f'Pairwise file list\n{PwFlLs}')
"""Create dataframe where focal species gene is has pair in all other species
and another dataframe that contains genes missing in atleast one species.
"""
CommonDf = pd.DataFrame()
for FlName in PwFlLs:
    #check if file name contains species identifiers at the right positions
    PwSpLs = FlName.split("_")[:2]
    SpFound = 0
    if PwSpLs[0] == FocalId: 
        OtherSpId = PwSpLs[1]
        if OtherSpId in SpIdDf.id.values: SpFound = 1
    elif PwSpLs[1] == FocalId: 
        OtherSpId = PwSpLs[0]
        if OtherSpId in SpIdDf.id.values: SpFound = 1
    if SpFound == 0: sys.exit("""\nOther species not at the first or second postion in the file name:\n
                              {}\n""".format(FlName))
    
    #merge pairwise dataframes to get a common dataframe
    PwDf = DropDuplicates("{}/{}".format(SpSynDir, FlName), FocalId, OtherSpId)
    PwDf.to_csv("{}/DupDropd_{}".format(OutFocalSpDataDir, FlName), sep='\t', header=True, index=True)
    CommonDf[OtherSpId] = PwDf.loc[:,OtherSpId].copy(deep=True)    
    
CommonIndex = CommonDf.dropna(axis=0, how='any').index    
NotCommonDf =  CommonDf.loc[CommonDf.index.difference(CommonIndex),:].copy(deep=True)
CommonDf.dropna(axis=0, how='any', inplace=True)
    
#create MSA
start_time = time.time()
MSA(CommonDf)
print("\n--- %s seconds ---" % (time.time() - start_time))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    