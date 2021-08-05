#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step - 4 Synteny and pairwise selection analysis

>>>> Check if the given species have same two letters identifiers,
    if so then the spcies identifier in the intermediate files will be more than
    two uppercase letters.

Given a list of species,
    conduct synteny assesment,
    run paml for each collinear gene pair for each species pair
    count sequence length of both proteins, 
    from the peptide alignment count, number of identical residues and substitutions. 
"""

import os, sys, pandas as pd, numpy as np, scipy
from Bio import AlignIO
import argparse
from Bio import SeqIO
import statsmodels.stats.multitest as multi
import scipy.stats
#from Bio import codonalign

######################  get the input species names and directory #######################


parser = argparse.ArgumentParser(description='Synteny compuration and pairwise selection analysis', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-s', '--species', required = True, help = 'Name of the species, without space and seperated by comma, for which pairwise blastp will be done')
parser.add_argument('-b', '--blastp_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/blastp', help = 'folder that has blastp files')
parser.add_argument('-g', '--gff_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/edited_gff', help = 'folder containing edited gff files')
parser.add_argument('-f', '--fasta_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/primTran', help = 'folder containing fasta files')
parser.add_argument('-B', '--blastp_suf', default = 'blastp.txt')
parser.add_argument('-T', '--trans_suf', default = 'transcript.fa')
parser.add_argument('-P', '--pep_suf', default = 'peptide.fa')
parser.add_argument('-G', '--gff_suf', default = 'gff3')
parser.add_argument('-d', '--doc_fol', default = '/groups/prodiv/projects/2019-10-21_mammals/doc/synteny', help = 'folder that will have the job submission files')
parser.add_argument('-p', '--paml_doc', default = '/groups/prodiv/projects/2019-10-21_mammals/doc/PairwisePaml', help = 'has control files') 
parser.add_argument('-o', '--out_fol', default = '/groups/prodiv/projects/2019-10-21_mammals/data/synteny', help = 'folder that will have the output')
parser.add_argument('-S', '--selection', default = 0)
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

def FlCheck(Dr, Prefix, Suffix, ContinueIfFileExists=False):
    """Check if file exists, return a string"""
    DrCheck(Dr, make=False)#check if directory exists
    FlLs = os.popen(r"""ls {} | grep '{}.*{}$'""".format(Dr, Prefix, Suffix)).read().strip('\n').split('\n')#grep file that ends with the suffix
    FlLs = list(filter(None, FlLs))#check if file exists, filter allows making list of length zero
    if len(FlLs) == 1:# file not in directory 
        if ContinueIfFileExists == True: return('continue')
        else: return(Dr+"/"+FlLs[0])#return file name
    elif len(FlLs) == 0:# file not in directory 
        sys.exit("""\n File not found.\n
                     Directory : {} \n
                     Prefix: {} \n
                     Suffix: {} \n""".format(Dr, Prefix, Suffix))
    elif len(FlLs) > 1:#multiple files, exit
        sys.exit("\n File identifier not unique.\nFiles; {}".format(FlLs))

def ParseColFl(ColinFl, TandemDupPairFl):
    """ 1. Sort collinear segments according to their score, rank them.
    Pass the rank to gene pairs. 
    2. Remove tandem duplicates from further analysis
    """
    #sort
    ColSorLs = os.popen(r""" grep '^## Alignment' {} |
                        sed 's/: score=/ /g' | sort -h -r -k4,4 |
                        sed 's/## Alignment //g' | sed 's/e_value=//g' |
                        sed 's/N=//g' """.format(ColinFl)).read().strip().split('\n') 
    #print(ColinFl)
    ColSorDf = pd.DataFrame([x.split() for x in ColSorLs], 
                             columns = ['alignment', 'score', 'e_value', 'N', 'chr', 'orient'])
    ColSorDf.alignment = ColSorDf.alignment.astype(int)
    ColSorDf.N = ColSorDf.N.astype(int)
    ColSorDf.score = ColSorDf.score.astype(float)
    ColSorDf.e_value = ColSorDf.e_value.astype(float)
    #Find species order in the collinearity file
    Id1 = ColSorDf.chr[0].split('&')[0][0:SpStrLetters+1]
    Id2 = ColSorDf.chr[0].split('&')[1][0:SpStrLetters+1]
    #load collinealirty file
    ColDf = pd.read_csv(ColinFl, comment = '#', sep = '\t', header = None,
                        names = ['No', Id1, Id2, 'e_value'])
    #rank gene pairs according to collinearity score
    ColDf['align'] = ColDf['No'].apply(lambda x: int(x.split('-')[0].strip()))
    RankSer = np.unique(ColDf['align'].values)# series with unique alignments values
        #From ColSorDf get the index number of every alignment and add one to make it the rank
    RankSer = pd.Series([ColSorDf.index[ColSorDf.alignment == x][0] + 1 for x in RankSer], index=RankSer) 
        #ColDf['rank'] = ColDf['align'].apply(lambda x: ColSorDf.index[ColSorDf.alignment == x][0] + 1) <<<< This takes time
    ColDf['rank'] = ColDf['align'].apply(lambda x: RankSer[x])
    ColDf.drop('align', axis=1, inplace=True)
    
    #drop Tandem Duplicats
    TanDupPairDf = pd.read_csv(TandemDupPairFl, sep=',', header=None)
    TanDupSer = TanDupPairDf.stack().reset_index(drop=True).drop_duplicates()
    ColDf = ColDf[~ColDf[Id1].isin(TanDupSer.values)]
    ColDf = ColDf[~ColDf[Id2].isin(TanDupSer.values)]
    return(ColDf)
    
    

def GetBestBlastMatch(BlastFl, gene):
    """Best reiprocal match based on match score
    """
    BestMatch = os.popen(r"""grep {} {} | sort -hr -k12,12 | head -n 1 | cut -f2
             """.format(gene, BlastFl)).read().strip('\n')
    return(BestMatch)
def BestReciprocalBlastPair(sp1, sp2, gene1, gene2):
    """Find if gene1 and gene2 are best reciprocal blast pair or not
    """
    BlastpFl1 = FlCheck(Dr = BlastDr, Prefix = "{}_vs_{}".format(sp1, sp2), Suffix = args["blastp_suf"])
    BlastpFl2 = FlCheck(Dr = BlastDr, Prefix = "{}_vs_{}".format(sp2, sp1), Suffix = args["blastp_suf"])
    if gene2 != GetBestBlastMatch(BlastpFl1, gene1): return('Not')
    if gene1 != GetBestBlastMatch(BlastpFl2, gene2): return('Not')
    else: return('Yes')
    
def PairwiseAnalysis(ColDf, sp1, sp2, sp1Id, sp2Id, PepFlSp1, PepFlSp2, TransFlSp1, TransFlSp2):
    """ Take gene names, transcript and, peptide file.
    Do selection analysis and count matches.
    """
    #copy control files
    os.popen(r"cp {}/*.ctl .".format(args["paml_doc"]))
    
    PepSp1 = SeqIO.index(PepFlSp1, 'fasta')
    PepSp2 = SeqIO.index(PepFlSp2, 'fasta')
    TransSp1 = SeqIO.index(TransFlSp1, 'fasta')
    TransSp2 = SeqIO.index(TransFlSp2, 'fasta')
    aaSp1 = '{}_aa'.format(sp1Id)
    aaSp2 = '{}_aa'.format(sp2Id)
    alSp1 = '{}_al'.format(sp1Id) 
    alSp2 = '{}_al'.format(sp2Id)
    AllAlignmntDf = pd.DataFrame()

    
    n = 0
    for indx, row in ColDf[[sp1Id, sp2Id]].iterrows():#iterate through the rows of the ColinDf to extract gene pairs
        #write tree file
        with open('Tree.txt', 'w') as fh: fh.write("({}, {});".format(row[sp1Id], row[sp2Id]))
        #create pep and trancript fasta files for the gene pair
        with open('protein1.fasta', 'w') as fh: 
            PepSeq1 = PepSp1[row[sp1Id]].seq
            SeqIO.write(PepSp1[row[sp1Id]], fh, 'fasta')
            if PepSeq1[-1] == '*': ColDf.at[indx,aaSp1] = len(PepSeq1) -1
            else: ColDf.at[indx,aaSp1] = len(PepSeq1)
        with open('protein2.fasta', 'w') as fh: 
            PepSeq2 = PepSp2[row[sp2Id]].seq
            SeqIO.write(PepSp2[row[sp2Id]], fh, 'fasta')
            if PepSeq2[-1] == '*': ColDf.at[indx,aaSp2] = len(PepSeq2) -1
            else: ColDf.at[indx,aaSp2] = len(PepSeq2)
        with open('transcripts.fasta', 'w') as fh: 
            SeqIO.write(TransSp1[row[sp1Id]], fh, 'fasta')
            SeqIO.write(TransSp2[row[sp2Id]], fh, 'fasta')
                
        #run strechter on the peptide file, fix gap open and gap extend
        os.popen(r"""    
        stretcher -asequence protein1.fasta -bsequence protein2.fasta -outfile Aligned_proteins.fasta -aformat fasta -gapopen 12 -gapextend 1
        """).read()
        
        #alignment parsing
        alignmnt = AlignIO.read('Aligned_proteins.fasta', 'fasta')
        alignDf = pd.DataFrame(np.array([list(rec) for rec in alignmnt])).transpose()
        overlapDf = alignDf[~alignDf.eq('-').any(1)].copy(deep=True)#drop rows that have gaps
        identityDf = overlapDf[overlapDf[0] == overlapDf[1]].copy(deep=True)
        IdRes = len(identityDf)
        Overlap = len(overlapDf)
        ColDf.at[indx,'Overlap'] = Overlap
        ColDf.at[indx,'IdRes'] = IdRes
        ColDf.at[indx,'Subs'] = Overlap - IdRes
        
        
        # percentage identity check and dropping spurious matches
        ColDf.at[indx,'Overlap_Id'.format(sp1Id)] = (IdRes*100.00)/Overlap
        Sp1Identity = (IdRes*100.00)/ColDf.at[indx,aaSp1]
        Sp2Identity = (IdRes*100.00)/ColDf.at[indx,aaSp2]
        ColDf.at[indx,'{}_Id'.format(sp1Id)] = Sp1Identity
        ColDf.at[indx,'{}_Id'.format(sp2Id)] = Sp2Identity
        if Sp1Identity < 95 or Sp2Identity < 95: #if not the best match drop the column
            if BestReciprocalBlastPair(sp1, sp2, ColDf.loc[indx, sp1Id], ColDf.loc[indx, sp2Id]) == 'Not':
                ColDf.drop(indx, axis=0, inplace=True)
                continue
                
        
        #alignment writing
        AlnRec = SeqIO.index('Aligned_proteins.fasta', 'fasta')
        AllAlignmntDf.loc[indx, sp1Id] = row[sp1Id]
        AllAlignmntDf.loc[indx, sp2Id] = row[sp2Id]
        AllAlignmntDf.loc[indx, alSp1] = str(AlnRec[row[sp1Id]].seq)
        AllAlignmntDf.loc[indx, alSp2] = str(AlnRec[row[sp2Id]].seq)
        
        if int(args["selection"]) == 1:
            os.popen(r"""
            pal2nal.pl Aligned_proteins.fasta transcripts.fasta -output paml -nogap> test.codon
            codeml codeml_for_neutral_selection.ctl
            codeml codeml_one_omega_for_all_branches.ctl
            """).read()
    
            #Parse PAML
            PamlOneOmegaFile = "result_same_omega_for_all_branches.txt"
            PamlNeutralOmegaFile = "result_neutral_selection.txt" 
            dN = os.popen(r"grep 'tree length for dN:' %s"%PamlOneOmegaFile).read().rstrip('\n').split()
            dS = os.popen(r"grep 'tree length for dS:' %s"%PamlOneOmegaFile).read().rstrip('\n').split()
            omega = os.popen(r"grep 'omega (dN/dS) =' %s"%PamlOneOmegaFile).read().rstrip('\n').split()
            lnL = os.popen(r"grep 'lnL(ntime:' %s"%PamlOneOmegaFile).read().rstrip('\n').split()
            lnL_Neutral = os.popen(r"grep 'lnL(ntime:' %s"%PamlNeutralOmegaFile).read().rstrip('\n').split()
            NS = os.popen(r"grep '^   3\.\.1' %s"%PamlOneOmegaFile).read().split('\n')[-2].split()
            if len(dS) == 0 or len(lnL_Neutral) == 0: 
                # In case PAML result file does not have the full results - skip
                ColDf.at[indx,'dS'] = np.NaN
                continue
            ColDf.at[indx,'dN'] = float(dN[-1])
            ColDf.at[indx,'dS'] = float(dS[-1])
            ColDf.at[indx,'NS'] = round(float(NS[2]))
            ColDf.at[indx,'SS'] = round(float(NS[3]))
            ColDf.at[indx,'Omega'] = float(omega[-1])
            ColDf.at[indx,'lnL'] = float(lnL[-2])
            ColDf.at[indx,'lnL_Neutral'] = float(lnL_Neutral[-2])
            #calculate LRT    
            LRT = 2*(float(lnL[-2]) - float(lnL_Neutral[-2]))
            P_Value = scipy.stats.distributions.chi2.sf(LRT,1)
            ColDf.at[indx,'P_Value'] = P_Value
#        print(n)
        # if n > 100: break
        # n += 1
    
    #do fdr correction    https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
    if int(args["selection"]) == 1: ColDf.loc[ColDf.P_Value.dropna().index, 'FdrAdjPvalue'] = multi.multipletests(ColDf.P_Value.dropna(), method = 'fdr_bh')[1]
    ColDf.to_csv('{}/{}_{}_CollinearDf.tsv'.format(OutSynDataDir, sp1Id, sp2Id), sep='\t', index=False, header=True)
    AllAlignmntDf.to_csv('{}/{}_{}/{}_{}.aln.tsv'.format(OutSynDataDir, sp1, sp2, sp1Id, sp2Id), sep='\t', index=False, header=True)
    
def SpeciesIdentifier(SpList, SpStrLetters = 1):
    """Take a list of species, pass a series with species names and unique
    Ids as the identifiers.
    """
    SpIdSer = pd.Series()
    for sp in SpList:
        sp = sp.strip() #remove leading and trailing spaces
        spId = sp.split('_')[0][0].upper() + sp.split('_')[-1][0:SpStrLetters].upper()
        if spId in SpIdSer.index:#if sp id is not unique includ more letters in the identifier
            SpStrLetters += 1
            SpIdSer, SpStrLetters = SpeciesIdentifier(SpList, SpStrLetters)
            return(SpIdSer, SpStrLetters)
        else: SpIdSer[spId] = sp
    SpIdSer.to_csv('{}/Id.txt'.format(OutSynDataDir), header = False, sep = '\t')
    return(SpIdSer, SpStrLetters)

"""Input"""
SpList = args["species"].split(',')
if len(SpList) < 2: sys.exit("""\n Less than two species given as input """)

FastaDr = args["fasta_dir"]
DrCheck(FastaDr)
BlastDr = args["blastp_dir"]
DrCheck(BlastDr)
GffDr = args["gff_dir"]
DrCheck(GffDr)
DocDr = args["doc_fol"]
DrCheck(DocDr)
OutDataDr = args["out_fol"]
DrCheck(OutDataDr)



"""Output"""
AllSpStr = args["species"].replace(' ', '_')
OutSynDataDir = OutDataDr + "/" + AllSpStr
DrCheck(OutSynDataDir, make = True)
OutSynDocDir = DocDr + "/" + AllSpStr
DrCheck(OutSynDocDir, make = True)


#make species identifier series and grab the number of species letters required to make unique ID
SpIdSer, SpStrLetters = SpeciesIdentifier(SpList)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>write species id ser to a file



"""Iterate over the species list"""
Run=0
Sp2List = SpList
for sp1 in SpList:
    sp1 = sp1.strip() #remove leading and trailing spaces
    sp1Id = SpIdSer[SpIdSer == sp1].index[0]
    sys.stdout.write("""\n Species 1:           '{}'
                     \n species identifier:     '{}'\n""".format(sp1, sp1Id))
    PepFlSp1 = FlCheck(Dr = FastaDr, Prefix = sp1, Suffix = args["pep_suf"])
    TransFlSp1 = FlCheck(Dr = FastaDr, Prefix = sp1, Suffix = args["trans_suf"])
    GffFlSp1 = FlCheck(Dr = GffDr, Prefix = sp1, Suffix = args["gff_suf"])
    TempGffFlSp1 = OutSynDataDir + "/" + sp1Id + ".gff"
    if Run == 0: #make sp1 gff file
        os.system(r"grep $'\tmRNA\t' " + GffFlSp1 + r""" | 
                awk '{print ($1, $4, $5, $9)}' | 
                sed -E 's/ID=transcript:([a-zA-Z0-9]+);.*$/\1/g' | 
                awk '{OFS = "\t"} {print ($1, $4, $2, $3)}' """ + r""" | 
                sed 's/^/{}/g' | sort -u >> {}""".format(sp1Id, TempGffFlSp1))
    #self blast file
    Sp1SelfBlastFl = FlCheck(Dr = BlastDr, Prefix = "{}_vs_{}".format(sp1, sp1), Suffix = args["blastp_suf"])
    os.popen(r"cp {} {}/{}.blast".format(Sp1SelfBlastFl, OutSynDataDir, sp1Id))
    #drop the species from Sp2list
    Sp2List = [x for x in Sp2List if x is not sp1]
    for sp2 in Sp2List:
        sp2 = sp2.strip() #remove leading and trailing spaces
        sp2Id = SpIdSer[SpIdSer == sp2].index[0]
        sys.stdout.write("""\n Species 2:           '{}'
                     \n species identifier:     '{}'\n""".format(sp2, sp2Id))
        PepFlSp2 = FlCheck(Dr = FastaDr, Prefix = sp2, Suffix = args["pep_suf"])
        TransFlSp2 = FlCheck(Dr = FastaDr, Prefix = sp2, Suffix = args["trans_suf"])
        BlastpFl1 = FlCheck(Dr = BlastDr, Prefix = "{}_vs_{}".format(sp1, sp2), Suffix = args["blastp_suf"])
        sys.stdout.write("""\n Blastp file 1: '{}'""".format(BlastpFl1))
        BlastpFl2 = FlCheck(Dr = BlastDr, Prefix = "{}_vs_{}".format(sp2, sp1), Suffix = args["blastp_suf"])
        sys.stdout.write("""\n Blastp file 2: '{}'""".format(BlastpFl2))
        Sp1Sp2DataDr = OutSynDataDir + "/{}_{}".format(sp1, sp2)
        DrCheck(Sp1Sp2DataDr, make = True)
        #concatenate blastp files
        BlastFl = "{}/{}_{}.blast".format(Sp1Sp2DataDr, sp1Id, sp2Id)
        Sp2SelfBlastFl = FlCheck(Dr = BlastDr, Prefix = "{}_vs_{}".format(sp2, sp2), Suffix = args["blastp_suf"])
        os.popen(r"cat {} {} {} {} > {}".format(BlastpFl1, BlastpFl2, Sp1SelfBlastFl, Sp2SelfBlastFl, BlastFl)).read()
        GffFlSp2 = FlCheck(Dr = GffDr, Prefix = sp2, Suffix = args["gff_suf"])
        TempGffFlSp2 = OutSynDataDir + "/" + sp2Id + ".gff"
        if Run == 0: #make sp2 gff file
            os.system(r"grep $'\tmRNA\t' " + GffFlSp2 + r""" | 
                    awk '{print ($1, $4, $5, $9)}' | 
                    sed -E 's/ID=transcript:([a-zA-Z0-9]+);.*$/\1/g' | 
                    awk '{OFS = "\t"} {print ($1, $4, $2, $3)}' """ + r""" | 
                    sed 's/^/{}/g' | sort -u >> {}""".format(sp2Id, TempGffFlSp2))
        #concatenate gff files
        GffFl = "{}/{}_{}.gff".format(Sp1Sp2DataDr, sp1Id, sp2Id)
        os.popen(r"cat {} {} > {}".format(TempGffFlSp1, TempGffFlSp2, GffFl)).read()
        #run MCScanX collinearity
        os.popen(r"/groups/prodiv/software/MCScanX/MCScanX {}/{}_{} -b 2 -m 2 -s 3".format(Sp1Sp2DataDr, sp1Id, sp2Id)).read()
        ColinFl = "{}/{}_{}.collinearity".format(Sp1Sp2DataDr, sp1Id, sp2Id)
        #run MCScanX tandem duplication
        TandemDupFl = "{}/{}_{}.TanDup.txt".format(Sp1Sp2DataDr, sp1Id, sp2Id)
        os.popen(r"""/groups/prodiv/software/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g {} -b {} -c {} -o {}""".format(
                GffFl, BlastFl, ColinFl, TandemDupFl)).read()
        TandemDupPairFl="{}/{}_{}.tandem".format(Sp1Sp2DataDr, sp1Id, sp2Id)
        #call function for parsing the collinearity file
        ColDf = ParseColFl(ColinFl, TandemDupPairFl)
        #run pairwise selection
        PairwiseAnalysis(ColDf, sp1, sp2, sp1Id, sp2Id, PepFlSp1, PepFlSp2, TransFlSp1, TransFlSp2)
        
        #call function for extracting genes within syntenic region but without collinear gene
        

    Run += 1
#delete temp blast and gff files
os.popen(r"rm -rf {}/*.gff".format(OutSynDataDir)).read()
os.popen(r"rm -rf {}/*.blast".format(OutSynDataDir)).read() 
os.popen(r"rm -rf {}/*/*.gff".format(OutSynDataDir)).read()      
os.popen(r"rm -rf {}/*/*.blast".format(OutSynDataDir)).read()  
    
    



#
#collinearity html parsing
#table = pd.read_html('./HS_vs_GG.html/hs1.html')
#MyTable = pd.DataFrame(table[0]).drop(0, axis =0) #select the first table on the page, though the page has only one table.
#MyTable.dropna(axis=1, how='all', inplace=True)
#
