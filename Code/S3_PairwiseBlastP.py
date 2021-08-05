#! /usr/bin/env python3
"""
Step -3 pairwise blastp

Accept 4 species name from input,
check if pairwise blastp files already exist.
else, check if corresponding peptide fasta files exits
    if yes, make blast database if not made already
    create jobsubmission script, 
    to run pairwise blastp for all possible combinations that do not exits already
    
"""

import os, sys
import argparse
 
#######################  get the input species names and directory #######################


parser = argparse.ArgumentParser(description='Does pairwise blastp for given species ', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--in_fol', default = '/groups/prodiv/projects/2019-10-21_mammals/data/primTran', help = 'folder that has original cds or petide fasta file')
parser.add_argument('-s', '--species', required = True, help = 'Name of the species, without space and seperated by comma, for which pairwise blastp will be done')
parser.add_argument('-o', '--out_fol', default = '/groups/prodiv/projects/2019-10-21_mammals/data/blastp', help = 'folder that has original cds or petide fasta file')
parser.add_argument('-F', '--pep_suf', default = 'peptide.fa')
parser.add_argument('-b', '--blastp_suf', default = 'blastp.txt')
parser.add_argument('-D', '--makeblast_suf', default = 'peptide.fa.phd')
parser.add_argument('-d', '--doc_fol', default = '/groups/prodiv/projects/2019-10-21_mammals/doc/blastp', help = 'folder that will have the job submission files')

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

def FlCheck(Dr, Prefix, Suffix, ContinueIfFileExists=False):
    """Check if file exists, return a string"""
    DrCheck(Dr, make=False)#check if directory exists
    FlLs = os.popen(r"""ls {} | grep '{}.*{}$'""".format(Dr, Prefix, Suffix)).read().strip('\n').split('\n')#grep file that ends with the suffix
    FlLs = list(filter(None, FlLs))#check if file exists
    if len(FlLs) == 1:#return file name
        if ContinueIfFileExists == True: return('continue')
        else: return(Dr+"/"+FlLs[0])
    elif len(FlLs) == 0:# file not in directory, 
        return('FileMissing')
    elif len(FlLs) > 1:#multiple files, exit
        sys.exit("\n File identifier not unique.\nFiles; {}".format(FlLs))


SbatchScript ="""#!/bin/bash
#  specify the job name
#SBATCH --job-name=XXX
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=72:00:00
#  maximum requested memory
#SBATCH --mem=4G
#  write std out and std error to these files
#SBATCH --error=job.%j.err
#SBATCH --output=job.%j.out
#  send a mail for job start, end, fail, etc.di
#SBATCH --mail-type=ALL
#SBATCH --mail-user=prabh@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard
# module load java/x64/8u121
"""

"""Input"""
SpList = args["species"].split(',')
if len(SpList) < 2: sys.exit("""\n Less than two species given as input """)
InDr = args["in_fol"]
DrCheck(InDr)
OutDr = args["out_fol"]
DocDr = args["doc_fol"]
DrCheck(OutDr)


"""Iterate over the species list"""
Run=0
for sp in SpList:
    sp = sp.strip() #remove leading and trailing spaces
    sys.stdout.write("""\n Species 1: '{}'""".format(sp))
    PepFl = FlCheck(Dr = InDr, Prefix = sp, Suffix = args["pep_suf"])
    sys.stdout.write("""\n Peptide fasta file for species 1 is: {} \n""".format(PepFl))
    MakeBlastDbCommand = ""
    if Run == 0:
        #run makeblastdb command if database is missing
        Check = FlCheck(Dr = InDr, Prefix = sp, Suffix = args["makeblast_suf"])
        if Check == 'FileMissing': os.popen(r"makeblastdb -in {} -dbtype prot -parse_seqids -hash_index".format(PepFl)).read()
    for sp2 in SpList:
        sp2 = sp2.strip() #remove leading and trailing spaces
        sys.stdout.write("""\n Species 2: '{}'""".format(sp2))
        PepFl2 = FlCheck(Dr = InDr, Prefix = sp2, Suffix = args["pep_suf"])
        sys.stdout.write("""\n Peptide fasta file for species 2 is: {} \n""".format(PepFl2))
        if Run == 0 and sp != sp2:
            Check = FlCheck(Dr = InDr, Prefix = sp2, Suffix = args["makeblast_suf"])
            if Check == 'FileMissing': os.popen("makeblastdb -in {} -dbtype prot -parse_seqids -hash_index".format(PepFl2)).read() 
        #check if blastp file already exists, if so continue
        Check = FlCheck(Dr = InDr, Prefix = sp2, Suffix = args["blastp_suf"])
        if Check == 'continue' : continue
        BlastpFaFl = OutDr + '/' + sp + '_vs_' + sp2 + '.' + args["blastp_suf"]
        BlastCommand = "blastall -i {}  -d {} -p blastp -e 1e-3 -b 5 -v 5 -m 8 -o {}".format(PepFl, PepFl2, BlastpFaFl)
        SbatchFl = "{}/Sbatch_{}_vs_{}_{}.sh".format(args["doc_fol"],sp, sp2, args["blastp_suf"])
        with open(SbatchFl, 'w+') as fh:
            fh.write(SbatchScript.replace("XXX", "{}_{}_blastp".format(sp, sp2)))
            fh.write(BlastCommand)
        os.system('sbatch {}'.format(SbatchFl))
    Run += 1
    

        
