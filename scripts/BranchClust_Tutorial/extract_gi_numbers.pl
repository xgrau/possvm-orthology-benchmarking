# ************************************************************
# Creation of Taxa Identification Table
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
# Changes added on June 25, 2007
# Original version assumed that number of digits in gi-number is equal to 8
# Now the number of digits in gi-number is estimated by the script. 
#

#!/usr/bin/perl -w

# check the location of fasta directory
# it should be in the same folder as the script

# 
# total number of digits in gi is estimated
#

open(LOG,"> gi_numbers.log") or die "can't open file gi_numbers.log\n";

sub assign_dots
{
#$TOTAL_DIGITS, $digits and $dots are global variables
($TOTAL_DIGITS, $digits)=@_;

my $k=$TOTAL_DIGITS-$digits;
$dots="";
while($k>0){
 $dots=$dots.".";
 $k--;
}
}

@files=();
@gis_table=();
@species_list=();

while(defined($file=glob("fasta/*.faa"))){

push(@files,$file);
@gis=(); 
 open(IN,"<$file") or die "can't open $file: $!\n";
   @lines=<IN>;
   $lines[0]=~/\[(.+?)\]/s;
   $species=$1;
   push @species_list, $species;
  
   for(my $k=0;$k<@lines;$k++){
   $line=$lines[$k];
   $gi_number="";
    if ($line=~/>gi/){
     $line=~/>gi\|(.+?)\|/s;
       $gi_number=$1;
       $nu_of_digits=length($gi_number);
       #initial precision is 4
       assign_dots($nu_of_digits,4);
       $FLAG_FOUND=0;
       for($j=0;$j<@gis;$j++){
        if($gi_number=~/$gis[$j]/){
          $FLAG_FOUND=1;
        }
       }
        if(!$FLAG_FOUND){
           @temp=split(//,$gi_number);
           $gi="";
           for($i=0;$i<$digits;$i++){
               $gi=$gi.$temp[$i];
           }
           $gi_template=$gi.$dots;
           push(@gis,$gi_template);
        }
    }
  }
  push @gis_table, [@gis];

close(IN);
}

close(OUT);

#--------------------------------------------------------
sub print_table
{  
  for($i=0;$i<$#gis_table+1;$i++){
   print LOG $species_list[$i]," |";
    for($j=0;$j<$#{$gis_table[$i]}+1;$j++){
      print LOG "\t$gis_table[$i][$j]";
    }
    print LOG "\n";
  }
}
#---------------------------------------------------------

#--------------------------------------------------------
sub print_table_to_file
{              
  for($i=0;$i<$#gis_table+1;$i++){
   print OUT $species_list[$i]," |";
    for($j=0;$j<$#{$gis_table[$i]}+1;$j++){
      print OUT "\t$gis_table[$i][$j]";
    }
    print OUT "\n";
  }
}
#---------------------------------------------------------


# we have now matrix with gi templates in @gis_table

sub find_doubles
{
my $i_max=$#gis_table+1;
 for($i=0;$i<$i_max;$i++){
    $j_max=$#{$gis_table[$i]}+1;
   for($j=0;$j<$j_max;$j++){
      $element=$gis_table[$i][$j];
      if(!defined($element)){
            die "$i,$j - element doesn't exist\n";
      }
        @indexes=();
        push @indexes,($i,$j);

          for($k=0;$k<$i_max;$k++){
            $m_max=$#{$gis_table[$k]}+1;
            for($m=0;$m<$m_max;$m++){            
              if(($element eq $gis_table[$k][$m])&&($k!=$i))
              { 
                print LOG "double entries: $element $species_list[$i] with $species_list[$k]\n";
                push @indexes, ($k,$m);
                
              }
           }
          }
      # if find doubles
      $n=@indexes;
      if($n > 2){
       return 1;
      }

   }
 }
@indexes=();
return 0;
}


sub refine
{
 my ($index1, $index2)=@_;
 my $gi_template=$gis_table[$index1][$index2];
 
 print LOG "REFINE: $gi_template of species $species_list[$index1] should be refined\n"; 

# calculate number of digits in the element
 $TOTAL_DIGITS=length($gi_template);
 my @temp=split(//,$gi_template);
 my $digits=$TOTAL_DIGITS;
 foreach (@temp){
  if($_ eq "\."){$digits--};
 }
 
 # if($digits==$TOTAL_DIGITS){die "two identical gi-numbers in the database:\n$gi_template for $species_list[$index1]\n"};

#increase number of digits to refine
 $digits++;

# delete this element from the table
 splice(@{$gis_table[$index1]},$index2,1);
 
 assign_dots($TOTAL_DIGITS,$digits);

 # open corresponding fasta file
 open(IN,"<$files[$index1]") or die "can't open file $files[$index1]\n";
 my @lines=<IN>; 
 close(IN);
 my @gis_refined=();
 for($k=0;$k<@lines;$k++){
  my $line=$lines[$k];
  my $gi_number="";
    if($line=~/>gi/){
      $line=~/>gi\|(.+?)\|/s;
      $gi_number=$1; 
       
        if($gi_number=~/$gi_template/){
          # extract first $digit numbers 
           @temp=split(//,$gi_number);
            $gi="";
            for($i=0;$i<$digits;$i++){
               $gi=$gi.$temp[$i];
            }
            # add it to the list if it is not yet included
            $gi_refined=$gi.$dots;
            my @found=grep(/$gi_refined/,@gis_refined);
            if(@found<1){
              push(@gis_refined,$gi_refined);
            }
        }
      
    }
 }
#update table
print LOG "REFINE: $gi_template is refined by ";
foreach(@gis_refined){
  print LOG "$_ ";
}
print LOG "\n";

push @{$gis_table[$index1]}, (@gis_refined);
}

while("doubles are found"){
   if(find_doubles()){
   
   for(my $i=0;$i<@indexes;$i=$i+2){
 
     $species1=$indexes[$i];
     $gene1=$indexes[$i+1];
   
     refine($species1,$gene1);
     print_table();
     print LOG "\n\n"

   }       
   }else{
    print LOG "No doubles found\n";
    open(OUT,">gi_numbers.out") or die "can't create file gi_numbers.out: $!\n";
    print_table_to_file;
    close(OUT);
    exit;
   }
}  
