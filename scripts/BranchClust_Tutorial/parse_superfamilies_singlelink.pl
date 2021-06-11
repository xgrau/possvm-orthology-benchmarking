# ************************************************************
# Assembles Superfamilies from Blast Results
# Single-linkage clustering
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
sub intersects
{
(my $line1, $line2)=@_;
 my @words1=split(/\s+/,$line1);
 my @words2=split(/\s+/,$line2);
 $common_word="";
 foreach $word1(@words1){
  if($word1 ne ""){
   my @comm=grep(/$word1/, @words2);
   if(@comm!=0) {
    $common_word=$word1;
    return 1;
   } 
  }
 }
 return 0;
}

# -------------------------------------------------------------------
sub combine_lines
{
(my $line1, $line2)=@_;
my @set1=split(/\s+/,$line1);
my @set2=split(/\s+/,$line2);
my @combined=();
push @combined, @set1;
  foreach $gi2 (@set2){
    if ($line1!~/$gi2/){
     push @combined, $gi2;
    }
  }
my $new_line=join(" ",@combined);
return $new_line;
}


# -------------------------------------------------------------------
# -------------------------------------------------------------------
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
 $logfile=$file;
 $outname=~s/.temp/.fam/;
 $logfile=~s/.temp/.log/;

 open(LOG,">$logfile") or die "can't open file $logfile\n";
 # print "processing $outname... \n";
 print LOG  "processing $outname... \n";
 
 open(IN,"<$file") or die "can't open $file: $! \n";
 open(OUT,">$outname") or die "can't open $outname: $! \n";
    @temp=<IN>;
  
 for($i=0;$i<@temp;$i++){
 # print "processing line $i\n";
   chomp($temp[$i]);
   $line1=$temp[$i];
   if($temp[$i] ne ""){
     print LOG "line $i combined with ";
     # print "line $i - \n";

     $found_intersect=1; 
     $common_word="";
     while($found_intersect)
     {     
       # print "cycling ...\n";
       $found_intersect=0; 
         for($j=0;$j<@temp;$j++){
           if($i!=$j){
             chomp($temp[$j]);
             $line2=$temp[$j];
             if($line2 ne ""){
                if(intersects($line1,$line2)){
                 #print "$i and $j intersects\n";
                 $new_line=combine_lines($line1,$line2);
                 $line1=$new_line;
                 $found_intersect=1; 
                 # delete the appended line from the array 
                 $temp[$j]="";
                 print LOG "$j($common_word) "; 
                 # print  "$j "; 
                } 
             }
           }
           }
      }#while
       
    #line $i is now consists of all intersected lines and can be printed to the file
    print OUT "$line1\n";   
    print LOG "\n";
    # print "\n";
    #line is cleaned from the array
    $temp[$i]=""; 
   }#if
 }#for i
 
#*****************************
 close(OUT);
 close(IN);
} #while

close(LOG);
