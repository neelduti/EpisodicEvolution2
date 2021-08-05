######## to extract alignment run the following command in terminal ######################

####initial alignment
grep Trancript_Id mafft.aln.tsv | tr '\t' '\n' > Output_Fast_File
#eg
grep ENSPTRT00000018053 mafft.aln.tsv | tr '\t' '\n' > ENSPTRT00000018053.mafft.fa

####Final filtered alignment

grep Trancript_Id gblocks.aln.tsv| tr '\t' '\n' > Output_Fast_File
#eg
grep ENSPTRT00000018053 gblocks.aln.tsv | tr '\t' '\n' > ENSPTRT00000018053.gb.fa
