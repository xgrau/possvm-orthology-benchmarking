# ************************************************************
# Blasts database composed of all genomes in question against itself
# 
# Part of BranchClust Tutorial
#
# BrunchClust algorithm
# version 1.01 November 10, 2006
# Copyright Maria Poptsova and J. Peter Gogarten 2006-2008
#
# E-mail: maria.poptsova@uconn.edu
# http://www.bioinformtics.org/branchclust/
#
# ************************************************************

#!/usr/bin/perl -w

#create 'blast' dir if it doesn't exist
if (!opendir(DIR,"blast")){
 mkdir("blast");
}else{
close(DIR);
}

$blast_input="fasta_all/allgenomes.faa";
$blast_output="blast/all_vs_all.out";

system("blastall -i $blast_input -d $blast_input -p blastp -o $blast_output -I T -e 1E-4 -F F -w 2 -a 2 -m 8");
    


