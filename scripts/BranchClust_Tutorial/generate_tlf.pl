# ************************************************************
# Generates TreeDyn Annotation File
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

if(!$ARGV[0]){die "usage: $0 CLUSTER_FILE.out.names>\n"};
$file=$ARGV[0];

$infile=$file;
$infile=~s/.names//;
open(IN1,"<$infile") or die "can't open file $infile\n";
$NU=0;
@cluster_nu=();
@family_nu=();
@lines=<IN1>;
for($i=0;$i<@lines;$i++){
  $line=$lines[$i];
  chomp($line);
  if($line=~/- CLUSTER -/){
   $NU++;
   $cluster=$lines[$i+1];
   @temp=split(/\s+/,$cluster);
   $cluster_nu[$NU]=scalar @temp;
   
  }
  if($line=~/- FAMILY -/){
    $family=$lines[$i+1];
    @temp=split(/\s+/,$family);
    $family_nu[$NU]=scalar @temp;
  }
}
close(IN1);

open(IN,"<$file") or die "can't open file $file\n";

$outfile="labelfile.tlf";
open(OUT,">$outfile") or die "can't open file $outfile\n";

sub search_table
{
 my $fam_gi=$_[0];
 for(my $i=0;$i<$#table+1;$i++){
  if($table[$i][0] eq $fam_gi){
   return $i;
 }
 }
return -1;
}

@table=();
$nu=0;
$flag_cluster=0;
$count=0;
while(defined($line=<IN>))
{
 chomp($line);
 if ($line=~/- CLUSTER -/){
  $nu++;
  $nu_in_cluster=0;
  $flag_cluster=1;
  next;
 }

if($line eq ''){
  $flag_cluster=0;
  next;
}

 if($flag_cluster){
   $nu_in_cluster++;
   # process line
   # each line consists of 3 part: gi, name, and species
   $line=~/gi\|(.+?)\|/s;
   $gi=$1;
   $line=~/\| (.+?) \[/s;
   $name=$1;
   $line=~/\[(.+?)\]/s;
   $taxa=$1;
   push @{$table[$count]},($gi,$nu,$name,$taxa);
   $count++;
 }

}

#refine information with families
seek(IN,0,0);
$flag_fam=0;
$flag_inpar=0;
$flag_outpar=0;
while(defined($line=<IN>)){
chomp($line);
 if ($line=~/FAMILY/){
   $flag_fam=1;
   next; 
 }
 if ($line=~/IN-PARALOGS/){
   $flag_inpar=1;
   next;
 }
 if ($line=~/OUT-OF_CLUSTER PARALOGS/){
   $flag_outpar=1;
   next; 
 }
 if($line eq ''){
  $flag_fam=0;
  $flag_inpar=0;
  $flag_outpar=0;
  next;
 }
 if(($flag_fam)||($flag_inpar)||($flag_outpar)){
   $line=~/gi\|(.+?)\|/s;
   $gi=$1;
   $i=search_table($gi);
   if($i==-1){die "error in family-match: gi=$gi\n"
   }else{
     if($flag_fam){
       push @{$table[$i]},"fam"; 
     }
     if($flag_inpar){
       push @{$table[$i]},"inpar";
     }
     if($flag_outpar){
       push @{$table[$i]},"outpar";
     }
   }
  }

  

} 
 

for($i=0;$i<$#table+1;$i++){
 $Nu_cluster=$table[$i][1];
 $nu_in_cluster=$cluster_nu[$Nu_cluster]; 
 $nu_in_fam=$family_nu[$Nu_cluster]; 
 print OUT "$table[$i][0] cluster $Nu_cluster family $table[$i][4] nu_cl $nu_in_cluster nu_fam $nu_in_fam name {$table[$i][2]} taxa {$table[$i][3]}\n";
 
}

close(IN);
close(OUT);


