# ************************************************************
# Creates one file from multiple fasta-files
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


#!usr/bin/perl -w

#create dir 
if (!opendir(DIR,"")){
 mkdir("fasta_all");
}else{
close(DIR);
}

system(" > fasta_all/allgenomes.faa");
while(defined($file=glob("fasta/*.faa")))
{
system("cat $file >> fasta_all/allgenomes.faa");
}
