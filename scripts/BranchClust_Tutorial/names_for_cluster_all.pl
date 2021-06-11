# ************************************************************
# Wrapper for names_for_cluster.pl
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

while(defined($file=glob("clusters/*.out")))
{
 print "processing  $file\n";
 `perl names_for_cluster.pl $file`;
}

