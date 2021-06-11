# ************************************************************
# Prints Name for Each Gene in a Cluster
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

if(!$ARGV[0]){
 die "eneter the name of the file to process \n";
}
$file=$ARGV[0];

open(IN,"< $file")or die "can't open file clusters.out $!\n";
@lines=<IN>;
$outfile=$file.".names";
open(OUT,"> $outfile") or die "can't open file $outfile: $! \n";

foreach $line(@lines){
 $line=~s/^\s+//g;  
 chomp($line);
 if($line ne "")
 {
    if(($line=~/CLUSTER/)or($line=~/FAMILY/)or($line=~/PARALOGS/)or($line=~/COMPLETE/)or($line=~/INCOMPLETE/)
     or($line=~/Number/)or($line=~/Run/)){
    print OUT $line,"\n";
 }else{ 
    $line=~s/-//g;
    $line=~s/^\s+//;  
    $line=~s/\s+/\,/g;
    $line=~s/\,$//;
     $out=`fastacmd -d fasta_all/allgenomes.faa -s $line -p T|grep '>'`;
     print OUT $out,"\n";
 } 
}
}

close(IN);
close(OUT);
