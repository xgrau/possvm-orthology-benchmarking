# ************************************************************
# Prints Summary for BranchClust Results
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

# identifiy number for "complete" set of species

$done=0;
while(defined($file=glob("clusters/*.list")))
{
 open(IN,$file) or die("can't open $file: $!\n");
 while(<IN>)
 {
  if($_=~/^COMPLETE: (.+?)$/){
   $total=$1;
   $done=1; 
   last;
  }
 }
close(IN);
if($done){last;}
}

 $incomplete=0;
 $complete=0;

for(my $i=$total-1;$i>0;$i--){
   $uncm[$i]=0;
}

while(defined($file=glob("clusters/*.list")))
{
 open(IN,$file) or die("can't open $file: $!\n");
 while(<IN>)
 {
  if($_=~/INCOMPLETE: /){
   $incomplete++;
  }
  if($_=~/COMPLETE: $total/){
   $complete++;
  }
  for(my $i=$total-1;$i>0;$i--){
    if ($_=~/INCOMPLETE: $i\n/){
      $uncm[$i]++;
    }
  }
 } 
 close(IN);
}

open(OUT,">summary.out") or die "can't open file summary.info\n";
 print OUT "complete: $complete\n";
 print OUT "incomplete: $incomplete\n";
 print OUT "total: ", $complete+$incomplete,"\n";
 print OUT "------ details -------\n";
 for(my $i=$total-1;$i>0;$i--){
    print OUT "incomplete $i: $uncm[$i]\n";
 }
close(OUT);

