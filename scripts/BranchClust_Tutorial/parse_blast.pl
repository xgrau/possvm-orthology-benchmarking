#************************************************************
# Parse Blast results
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

use lib "/Users/mpoptsova/bioperl-1.5-my";
use Bio::SearchIO;

#create 'parsed' if it doesn't exist
if(!opendir(DIR,"parsed")){
 mkdir("parsed");
}else{
close(DIR);
}

$infile="blast/all_vs_all.out";
$outfile="parsed/all_vs_all.parsed";

open (OUT, ">$outfile") || die "Cannot open file $outfile $!\n";

my $in = new Bio::SearchIO(-format => 'blasttable',
                           -file => "$infile");

# Because it is blast of a database of N genomes against itself,
# first hit for each gene is the gene itself.
# That is why we assemble only hits.

while(my $result = $in->next_result){
  while($hit = $result->next_hit()){
   # take only first hsp for every hit, it has the best e-value
         
   # exctract gene number
   $hit->name()=~/\|(.+?)\|/s;  
   $gene=$1;    
   print OUT  "$gene\t";    
  }
  print OUT "\n";
}
close (OUT);
