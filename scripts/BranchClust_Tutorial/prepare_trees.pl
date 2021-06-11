# ************************************************************
# Prepares Trees to use with BranchClust
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

# this program prepare clustalw generated trees for further analyisis 
# it converts tree into one line

while(defined($file=glob("dist/*.ph")))
{
  open(PH,"<$file") or die "can't open $file: $!\n";
  $outfile=$file.".my";
  open(TREE,">$outfile") or die "can't open $outfile: $!\n";
   @ph_tree=<PH>;
    for($i=0;$i<@ph_tree;$i++){
     chomp($ph_tree[$i]);     
    }
  print TREE @ph_tree;
  close(TREE);  
  close(PH);
}


while(defined($file=glob("dist/*.ph.my")))
{
   open(IN,"<$file") or die "can't open file $file: $!\n";
   @tree=<IN>;
   close(IN);
   $outfile=$file;
   $outfile=~s/.ph.my/.tre/;
   open(OUT,">$outfile") or die "can't open file $outfile: $!\n";
   $tree[0]=~s/gi\|//g;
   $tree[0]=~s/\|ref\|(.+?):/:/g;
   $tree[0]=~s/\|sp\|(.+?):/:/g;
   $tree[0]=~s/\|gb\|(.+?):/:/g;
   print OUT $tree[0];
   close(OUT);
}

#create trees if it doesn't exist
if(!opendir(DIR,"trees")){
 mkdir("trees");
}else{
close(DIR);
}

system("mv dist/*.tre trees/");
