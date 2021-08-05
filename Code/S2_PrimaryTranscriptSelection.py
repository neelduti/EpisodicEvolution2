#! /usr/bin/env python3
"""
Step -2 Primary Transcript Selection

Open the gff file,
Select a gene, then check how many mRNAs it is the parent of,
if there are two or more mRNAs,iterate over them and select the
CDSs theY are parents of.
add the length of CDSs and select the longest mRNA

Extract the transcripts in a new transcipt fasta file

Translate the transcripts to create peptide files
"""

import os, sys, pandas as pd, numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
 
#######################  get the input species names and directory #######################


parser = argparse.ArgumentParser(description='Selects longest isoform pergene and creates both transcript and peptide fasta files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-C', '--in_cds_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/cds', help = 'folder that has original cds or petide fasta file')
parser.add_argument('-F', '--cds_suf', default = 'fa')
parser.add_argument('-B', '--biotype', default = 'biotype=protein_coding')
parser.add_argument('-I', '--in_gff_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/gff', help = 'folder that has original gff file')
parser.add_argument('-G', '--gff_suf', default = 'gff3')
parser.add_argument('-Z', '--zip_suf', default = '.gz')
parser.add_argument('-s', '--species', required = True, help = 'Name of the species')
parser.add_argument('-g', '--gff_out_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/edited_gff', help = 'folder that will have edited gff files')
parser.add_argument('-f', '--fasta_out_dir', default = '/groups/prodiv/projects/2019-10-21_mammals/data/primTran', help = 'folder that will fasta files')
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
            os.mkdir(Dr)

def FlCheck(Dr, Prefix, Suffix, ZipSuf = '.gz', make=False):
    """Check if unzipeed file exists, if yes return it
    otherwise find the zipped version, unzip it and
    then return the unzipped file"""
    DrCheck(Dr, make=False)#check if directory exists
    FlLs = os.popen(r"""ls {} | grep '{}.*{}$'""".format(Dr, Prefix, Suffix)).read().strip('\n').split('\n')
    FlLs = list(filter(None, FlLs))#check if file exists
    if len(FlLs) == 1:#return file name
        return(Dr+"/"+FlLs[0])
    elif len(FlLs) == 0:# file not in directory, check if zipped file exists
        ZpFlLs = os.popen(r"""ls {} | grep '{}.*{}{}$'
                          """.format(
                          Dr, Prefix, Suffix, ZipSuf)
                          ).read().strip('\n').split('\n')#check if zipped file exists
        ZpFlLs = list(filter(None, ZpFlLs))
        if len(ZpFlLs) == 1:#zipped file foound
            sys.stdout.write("""\n Unzipping file: '{}'
                             in the directory '{}'.""".format(ZpFlLs[0], Dr))
            os.popen(r"""gzip -d {}""".format(Dr+"/"+ZpFlLs[0])).read()
            return((Dr + "/" + ZpFlLs[0].rstrip(ZipSuf)))
        elif len(ZpFlLs) == 0:#no zipped file either
            sys.exit("""\n Bad file identifier: \n 
                     Directory : {} \n
                     Prefix: {} \n
                     Suffix: {} \n
                     Zip: {} \n
                     Neiter zipped or unzipped files found.
                     """.format(Dr, Prefix, Suffix, ZipSuf))
        elif len(ZpFlLs) > 1: #multiple zipped files found
            sys.exit("\n Zipped file identifier not unique.\nFiles; {} \n".format(ZpFlLs))
    elif len(FlLs) > 1:#multiple files, exit
        sys.exit("\n File identifier not unique.\nFiles; {}".format(FlLs))

def GffParse(sp, GffFl, CdsFl, OutGffFl, OutTranFaFl, OutTranLs):
    """Parse gff file, 1st extract genes assosiated with single transcripts,
    select protin coding genes, and 
    for genes associated with multiple transcripts pict the longest one 
    as the primary transcript.
        Using the transcript identifier, create the primary transcript fasta 
    file.
    """
    sys.stdout.write("""\n Parsing Gff file{} """.format(GffFl))
    GffDf = pd.read_csv(GffFl, header=None, sep='\t', comment='#', low_memory=False)#read gff file as pandas dataframe
    GffDf = GffDf[(GffDf[2]=='mRNA') | (GffDf[2]=='gene') | (GffDf[2]=='CDS')]#reomve lines other than genes, mrna and CDS
    #create an identifiler column that has either the gene or the transcript identifier
    GffDf['Identifier'] = GffDf[~(GffDf[2] == 'CDS')][8].apply(lambda x: x.split('ID=')[1].split(';')[0].split(':')[1])
    GffDf.loc[GffDf[(GffDf[2] == 'CDS')].index,'Identifier'] = GffDf[(GffDf[2] == 'CDS')][8].apply(lambda x: x.split('Parent=transcript:')[1].split(';')[0]) 
    MrnaDf = GffDf[GffDf[2]=='mRNA'].copy(deep=True)#select the mRNA rows
    MrnaDf['gene'] = MrnaDf[8].apply(lambda x: x.split('Parent=gene:')[1].split(';')[0])#introduce gene column
    MrnaDf['tran'] = MrnaDf[8].apply(lambda x: x.split('ID=transcript:')[1].split(';')[0])#introduce transcript column
    #drop non-protein coding genes and mRNA
    # MrnaDf['biotype'] = MrnaDf[8].apply(lambda x: 'protein_coding' if args['biotype'] in x else 'other')#extract biotype
    MrnaDf['biotype'] = MrnaDf[8].apply(lambda x: x.split('biotype=')[1].split(';')[0])
    sys.stdout.write("""\nlength of mRNA Df : {}\nMrnaDf Head 1 \n{}\n""".format(len(MrnaDf), MrnaDf.head()))#
    MrnaDf = MrnaDf[(MrnaDf['biotype']=='protein_coding')]#drop non-protein coding mRNA
    GeneCnt = MrnaDf.gene.value_counts()#gene count
    sys.stdout.write("""\nlength of mRNA Df : {}\nMrnaDf Head 2 \n{}\n""".format(len(MrnaDf), MrnaDf.head()))
    GeneCnt_1 = GeneCnt[GeneCnt==1]# extract genes with single trancripts
    MrnaDf.set_index('gene', inplace=True, drop=False)#index by gene
    GeneDf = GffDf[GffDf[2]=='gene'].copy(deep=True)#select the gene rows
    GeneDf['biotype'] = GeneDf[8].apply(lambda x: x.split('biotype=')[1].split(';')[0])#extract biotype
    sys.stdout.write("""\nGeneDf Head 1 \n{}\n""".format(GeneDf.head()))#
    GeneDf = GeneDf[(GeneDf['biotype']=='protein_coding')]#drop non-protein coding genes
    sys.stdout.write("""\nGeneDf not protein Head 1 \n{}\n""".format(GeneDf[~(GeneDf['biotype']=='protein_coding')].head()))
    sys.stdout.write("""\nGeneDf Head 2 \n{}\n""".format(GeneDf.head()))#
    GeneDf['gene'] = GeneDf[8].apply(lambda x: x.split(';')[0].lstrip('ID=gene:'))#extract gene name
    GeneDf.set_index('gene', inplace=True, drop=False)#index by gene
    MrnaDf = MrnaDf.loc[MrnaDf.index.intersection(GeneDf.index).unique(),:]
    #select one out of the multiple mRNAs
    PrimMrnaDf = MrnaDf.loc[MrnaDf.index.intersection(GeneCnt_1.index),:].copy(deep=True)#
    MultiMrnaDf = MrnaDf.loc[MrnaDf.index.difference(GeneCnt_1.index),:].copy(deep=True)
   
    SelMrnaDf = pd.DataFrame()
    MaxTranSer = pd.Series(dtype=int)
    for Gene in MultiMrnaDf.index.unique():
        TranSer = pd.Series(dtype=int)
        for Tran in MultiMrnaDf.loc[Gene, 'tran']:
            tmpCDSDf = GffDf[(GffDf.Identifier == Tran) & (GffDf[2] == 'CDS')].copy(deep=True)
            CodLen = np.sum(tmpCDSDf[4] - tmpCDSDf[3])
            TranSer[Tran] = CodLen
        if TranSer.max() == 0: 
            sys.exit("""There is an error in selection of the longest transcript.
                     The maximum transcript lengh reads zero.
                     """)
        MaxTranSer[Tran] = TranSer.max()
        SelMrnaDf = SelMrnaDf.append(MultiMrnaDf[MultiMrnaDf['tran']==TranSer.idxmax()].copy(deep=True))
    PrimMrnaDf = PrimMrnaDf.append(SelMrnaDf)#Comple primary mRNA dataframe for all genes
    #create edited gff file and primary trancript list
    GffDf.set_index('Identifier', inplace=True, drop=False)
    PrimMrnaDf.set_index('Identifier', inplace=True, drop=False)
    GffDf = GffDf.loc[GffDf.index.intersection(PrimMrnaDf.index).unique(),:]#the gff file has same identifier s for CDSs and mRNA
    GffDf[8] = GffDf[8].str.replace('Parent=gene:', 'Gene=')
    pd.Series(PrimMrnaDf.index).to_csv(OutTranLs, header=False, index=False)#write gene name list
    #GffDf.to_csv(OutGffFl, header=False, index=False, sep='\t')
    #create edited fasta file and gff file
    os.popen(r"""grep -F -f {} {} > {}
             """.format(OutTranLs, GffFl, OutGffFl)).read()
    os.popen(f"sed -i 's/Parent=gene/Gene=/g' {OutGffFl}").read()
    sys.stdout.write("""\n Editing fasta file {}""".format(CdsFl))
    os.popen(r"""cp {} tmp.{}.cds.fa""".format(CdsFl, sp)).read()#copy cds fasta file
    os.popen(r"""sed -i -E 's/^(>[a-zA-Z0-9]+).[0-9]+/\1/g' """ + """tmp.{}.cds.fa""".format(sp)).read()
    os.popen(r"""perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV'""" + """ {} tmp.{}.cds.fa > {}""".format(
            OutTranLs, sp, OutTranFaFl)).read()
    os.popen(r"""rm tmp.{}* """.format(sp)).read()

def TransToPep(OutTranFaFl, OutPepFaFl):
    """Create peptide fasta file from transcript fasta file.
    """
    sys.stdout.write("""\n Translating fasta file {}'""".format(OutTranFaFl))
    records = SeqIO.index(OutTranFaFl, 'fasta')
    PepLs = []
    for Id in records.keys():
        PepLs.append(SeqRecord(records[Id].seq.translate(to_stop=False), id = Id))
        #PepLs.append(SeqRecord(records[Id].seq.translate(to_stop=True), id = Id))
    SeqIO.write(PepLs, OutPepFaFl, 'fasta')

"""Input"""
# SpList = args["species"].split(',')
sp = args["species"]
InCdsDr = args["in_cds_dir"]
DrCheck(InCdsDr)
InGffDr = args["in_gff_dir"]
DrCheck(InGffDr)
GffDr = args["gff_out_dir"]
DrCheck(GffDr)
FastaDr = args["fasta_out_dir"]
DrCheck(FastaDr)


"""Iterate over the species list or run for a single species"""

# for sp in SpList:
sys.stdout.write("""\n Species: '{}'""".format(sp))
sp = sp.strip() #remove leading and trailing spaces
CdsFl = FlCheck(Dr = InCdsDr, Prefix = sp, Suffix = args["cds_suf"])
sys.stdout.write("""\n Cds file is: {} \n""".format(CdsFl))
GffFl = FlCheck(Dr = InGffDr, Prefix = sp, Suffix = args["gff_suf"])
sys.stdout.write("""\n Gff file is: {} \n""".format(GffFl))
OutGffFl =    GffDr + '/' + sp + '.transcript.gff3'
OutTranFaFl = FastaDr + '/' + sp + '.transcript.fa'
OutTranLs = GffDr + '/' + sp + '.transcript.txt'
GffParse(sp, GffFl, CdsFl, OutGffFl, OutTranFaFl, OutTranLs)
OutPepFaFl = FastaDr + '/' + sp + '.peptide.fa'
TransToPep(OutTranFaFl, OutPepFaFl)    

