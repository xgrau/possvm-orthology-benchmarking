# ************************************************************
# Prints detailed summary from BranchClust results
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

#!/usr/bin/perl

# we will print text table in the format
# superfamily_##  family_##  nu_of_genes_in_the_family  nu_of_paralogs  family_name
#

open(OUT,">detailed_summary.out") or die "can't create file detailed_summary.out\n";
print OUT "superfamily_##\tfamily_##\tnu_of_genes_in_the_family\tnu_of_paralogs\tfamily_name\n";

open(FAM,">families-names.list") or die "can't create file: $!\n";

while(defined($file=glob("clusters/*.out")))
{
 $file=~/\/clusters_(.+?).out/s;
 $sfam=$1;
 $fam_count=0;
  open(IN,"<$file") or die "can't open file $file\n";
  @lines=<IN>;
    for(my $i=0;$i<@lines;$i++){
      if($lines[$i]=~/CLUSTER/){
        @cluster=();
        @cluster=split(/\s+/,$lines[$i+1]);
        $nu_in_cluster=scalar @cluster; 
      }
      if($lines[$i]=~/FAMILY/){
        @fam=();
        @fam=split(/\s+/,$lines[$i+1]);
        $nu_in_fam=0;
         foreach $member(@fam){
           if($member ne "-"){
              $nu_in_fam++;
              $gi=$member;
           }
         }
        $fam_count++;
         
    
        $nu_of_paralogs=0;
        $nu_of_paralogs=$nu_in_cluster-$nu_in_fam;
       
        $out=`fastacmd -d fasta_all/allgenomes.faa -s $gi -p T|grep '>'`;
        $out=~/\| (.+?) \[/s;
        $family_name=$1; 
        print OUT "$sfam\t$fam_count\t$nu_in_fam\t$nu_of_paralogs\t$family_name\n";
        if($nu_in_fam>=1){
          
          print FAM "$sfam.$fam_count.$nu_in_fam.$nu_of_paralogs\t@fam\n";   
        }
      }
    }
}

close(OUT);
close(FAM);
