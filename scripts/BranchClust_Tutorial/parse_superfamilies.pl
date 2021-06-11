# ************************************************************
# Assembles Superfamilies from Blast Results
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

open(LOG,">parse_superfam.log") or die "can't open file check.log \n";

if(!$ARGV[0]) {die "usage: $0 <nu_diff_gen>\n"};
#relax parse requirements that all genomes should be present
$nu_diff_gen=$ARGV[0];
# --------------------------------------------------------------------
# before running this script we must create a file gi_numbers.out with
# extract_gis.pl
 
 open(GIS,"< gi_numbers.out") or die "can't open file gi_numbers.out: $!\n";
 @gis_legend=<GIS>;
 close(GIS);

$nu_of_genomes=@gis_legend;
 
 for($i=0;$i<@gis_legend;$i++){
   chomp($gis_legend[$i]);
   # remove species name
   $gis_legend[$i]=~s/^(.+?)\|\t//;
   # print "$i $gis_legend[$i] \n";
   @words=();
   @words=split("\t",$gis_legend[$i]);
  #check for empty
   foreach(@words){
    if($_ eq "") {die "there are empty records in gi_numbers.out\n";}
   }
   push @gis_table, [@words];
 
}


# -------------------------------------------------------------------
sub search_intersections
{
  my ($current_gi, $current_index, $ref, $ref_combined)=@_;
 
   # COMBINE CLUSTERS WITH PARALOGS
   # search all occurence of the $current_gi
   # remember indexes of found clusters in all_clusters hash

   my %all_clusters=();
   
   for(my $j=0;$j<@$ref;$j++){
    chomp($$ref[$j]);
    if(($$ref[$j]=~/$current_gi/)&&($j!=$current_index)){
       my $len=length($$ref[$j]);
       print LOG "$current_gi found both in the line $current_index and $j\n";
       $all_clusters{$j}=$len;
     }
   }#for j

  
  my $size = scalar keys (%all_clusters);
  if ($size==0){
  return 0;
  }else{  

  
   foreach $index (keys %all_clusters){
     print LOG "$index ";
     # add index cluster to combined
     # check for duplicates
     my @check;
     chomp($$ref[$index]);
     my @genes=split("\t",$$ref[$index]);
        for(my $j=0;$j<@genes;$j++){
        # first check if an element is alredy included in the array:
        my @check=();
        my $check_gi=$genes[$j];
        @check = grep(/$check_gi/,@combined);
        # print "check: $genes[$j] @check\n";

        # if not found then add it to the combined cluster
        if (@check==0){
         push(@combined,$genes[$j]);
        }
    }

      $$ref[$index]="";
    }  #foreach index
 
 }

}
# -------------------------------------------------------------------
# convert parsed to temp
# each line of temp should contain more or equal to the defined number 
# >= $nu_diff_gen (parameter)
#

while(defined($file=glob("parsed/*.parsed"))){
 $outname=$file;
 $outname=~s/.parsed/.temp/;
 open(IN,"<$file") or die "can't open $file: $! \n";
 open(OUT,">$outname") or die "can't open $outname: $! \n";

    while(defined($line=<IN>)){
     chomp($line);

     $genome_count=0;
#--------------------------------------------------------------
      
    $i_max=$#gis_table+1;
    for($i=0;$i<$i_max;$i++){
      $j_max=$#{$gis_table[$i]}+1;
      for($j=0;$j<$j_max;$j++){
         if ($line=~$gis_table[$i][$j]){
          $genome_count++;
          last;
         }
      }
    }

#-------------------------------------------------------------      
     if ($genome_count>=$nu_diff_gen){
      print OUT "$line\n";
     }     
        
          
    }
 close(OUT);
 close(IN);
}

#-----------------------------------------------------------------------

#The file temp contains superfamily with number of different taxa equals to $nu_diff_gen 

#--------------------------------------------------------------------
while(defined($file=glob("parsed/*.temp"))){
 $outname=$file;
 $outname=~s/.temp/.fam/;
 print "processing $outname... \n";
 print LOG  "processing $outname... \n";
 
 open(IN,"<$file") or die "can't open $file: $! \n";
 open(OUT,">$outname") or die "can't open $outname: $! \n";
    @temp=<IN>;
    for($i=0;$i<@temp;$i++){
       chomp($temp[$i]);
       if ($temp[$i] ne ""){
         
         print LOG "cluster: $i ";

         @gis=split("\t",$temp[$i]);
            
           @combined=();
           push (@combined, @gis);

           foreach $current_gi (@gis)
           {
             search_intersections($current_gi,$i, \@temp,\@combined);
           
           }  
             
           #print out cluster;
           foreach $gene(@combined){
             print OUT "$gene\t";   
           }  
           print OUT "\n";
           print LOG "\n";

          $temp[$i]="";
  
       } #if
    } #for i
 close(OUT);
 close(IN);
} #while

close(LOG);
