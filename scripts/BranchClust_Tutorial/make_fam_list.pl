# ************************************************************
# Prints list of families up to the  <NU_OF_INCOMPLETE TAXA>
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

if((!$ARGV[0])or(!$ARGV[1])){die "usage: $0 <NU_OF_ALL_TAXA> <NU_OF_INCOMPLETE_TAXA>\n";}
$nu_compl=$ARGV[0];
$nu_incompl=$ARGV[1];

open(OUT,">families.list") or die "can't open file families.list\n";

while(defined($file=glob("clusters/*.list")))
{
 open(IN,$file) or die("can't open $file: $!\n");
 @lines=<IN>;
 for($i=0;$i<@lines;$i++)
 {
  $_=$lines[$i];
  if($_=~/COMPLETE: (.+?)$/){
    $nu=$1;
    if(($nu<=$nu_compl)&&($nu>=$nu_incompl)){
     print OUT "$lines[$i-1]";
    }
  }     
 } 
 close(IN);
}
close(OUT);
