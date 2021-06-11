# ************************************************************
# Prepares Fasta Files for Each Superfamily
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

if(!$ARGV[0]){die "usage: $0 <file-with-fam>"}
$families=$ARGV[0];

#we are currently in the /clustalw/directory
#create fa if it doesn't exist
if (!opendir(DIR,"fa")){
 mkdir("fa");
}else{
close(DIR);
}

#make sure you have database with all genomes in this location:
$all_database="fasta_all/allgenomes.faa";

 #open files with families;
  open(IN, "< $families") or die "cannot open $families: $!";
  $count=1;

  while(defined($line=<IN>))
  {
  $line=~s/\t/,/g;
  chomp($line);

 # we are in the clustalw directory
  $all_families=`fastacmd -d $all_database -s $line -p T`;
    
  $outfile='fa/fam_'.$count.'.fa';
  open(OUT, "> $outfile") or die "cannot open $outfile: $!";  

  print OUT $all_families;
  close OUT;
  $count++;
  }

  close(IN);


