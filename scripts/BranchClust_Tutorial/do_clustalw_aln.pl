# ************************************************************
# Superfamily alignment with clustalw
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

#! /usr/bin/perl -w

#create 'dist' if it doesn't exist
if (!opendir(DIR,"dist")){
 mkdir("dist");
}else{
close(DIR);
}


while(defined($filename=glob("fa/*.fa")))
{
  print "$filename\n";
  # clustalw each file
  system("clustalw -infile=$filename -align -type=protein");
}

system("mv fa/*.aln dist/$d");
 

