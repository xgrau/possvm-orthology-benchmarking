# ************************************************************ 
# BrunchClust algorithm wrapper
#
# Copyright Maria Poptsova and J. Peter Gogarten 2006-2008
#
# E-mail: maria.poptsova@uconn.edu
# http://www.bioinformtics.org/branchclust/ 
#
# ************************************************************
#!/usr/bin/perl

if(!$ARGV[0]){die "usage $0 <MANY>" }
$many=$ARGV[0];

if(!opendir(DIR,"clusters")){
 mkdir("clusters");
}else{
 close(DIR);
}

open(LOG,">branchclust.log") or die "can't open file branchclust.log\n";

while(defined($file=glob("trees/*.tre")))
{
 
 $file=~/trees\/fam_(.+?).tre/s;
 $nu=$1;
  $time=localtime();
  print LOG "tree $nu - $time\n";
 `perl branchclust.pl $file $many`;
 `mv clusters.out clusters/clusters_$nu.out`;
 `mv families.list clusters/families_$nu.list`;
 `mv cluster.log clusters/cluster_$nu.log`;
}
$time=localtime();
print LOG "finished: $time\n";
close(LOG);

