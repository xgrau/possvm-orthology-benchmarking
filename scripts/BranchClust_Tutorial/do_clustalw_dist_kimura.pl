# ************************************************************
# Tree Reconstruction
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

# Tree reconstruction, using distance method with kimura correction
# trees will be generated in the same directory 'dist' with extension *.ph

while(defined($filename=glob("dist/*.aln")))
{
 # clustalw each file
 system("clustalw -infile=$filename -tree -OUTPUTTREE=dist -kimura");
}
 
