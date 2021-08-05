#Run in terminal
ProFl=$(pwd)
################### 1. data download
mkdir $ProFl/data $ProFl/data/cds $ProFl/data/gff $ProFl/doc
cd $ProFl/data/cds
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-98/fasta/*/cds/*.cds.all.fa.gz .
cd $ProFl/data/gff
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-98/gff3/*/*.98.gff3.gz .
cd ..
gzip -d ./*/Homos_sapiens*
gzip -d ./*/Pan_troglodytes*
gzip -d ./*/Gorilla_gorilla*
gzip -d ./*/Nomascus_leucogenys*
#delete all the zipped files once unzipping is complete
rm ./*/*.gz

################### 2. Primary transcript selection 
#if possible run these as four separate jobs on a HPC
python3 $ProFl/code/S2_PrimaryTranscriptSelection.py -s Homos_sapiens
python3 $ProFl/code/S2_PrimaryTranscriptSelection.py -s Pan_troglodytes
python3 $ProFl/code/S2_PrimaryTranscriptSelection.py -s Gorilla_gorilla
python3 $ProFl/code/S2_PrimaryTranscriptSelection.py -s Nomascus_leucogenys

################### 3. Pairwise all vs all BLASTp
python3 $ProFl/code/S3_PairwiseBlastP.py -s Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys



################### 4 Pairwise synteny
mkdir $ProFl/data/synteny $ProFl/doc/synteny
#if possible run this job on a HPC
python3 $ProFl/code/S4_Synteny.py -s Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys


################### 5 MSA
mkdir $ProFl/data/msa $ProFl/doc/msa
python3 $ProFl/code/S5_MSA.py -s Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys -f Homo_sapiens


################### 6 Alignment parsing
python3 $ProFl/code/S6_ParseAlignment.py -f $ProFl/data/msa/Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys/Homo_sapiens/gblocks.aln.tsv


################### 7 RF Metric
cd $ProFl/data/msa/Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys/Homo_sapiens/Subs_gblocks.aln
mkdir WithNull
cd WithNull
python3 $ProFl/2019-10-21_mammals/code/S7a_RFmetricSE.py
python3 $ProFl/2019-10-21_mammals/code/S7b_Ztest.py -f PoissonCorrectedDf.tsv

################### 8 Simulation and candidate identification
cd $ProFl/data/msa/Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys/Homo_sapiens/Subs_gblocks.aln
python3 $ProFl/2019-10-21_mammals/code/S8a_Simulation.py
python3 $ProFl/2019-10-21_mammals/code/S8b_Candidates.py

################### 9 Figures and tables
cd $ProFl/data/msa/Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys/Homo_sapiens
mkdir Figures
cd Figures
python3 $ProFl/2019-10-21_mammals/code/Figures.py
