#! /usr/bin/perl -w
#

use strict;

my $file = shift @ARGV;
open my $fh, '<', $file or die "ERROR: Cannot open $file";

#total rattle, pt events, and transition events
my ($num_rat, $num_pt, $num_trans);

#index number from the last step
my ($prev_index, $prev_step);

#proton transfer
my $pt = 0;

while (my $line = <$fh>){
   chomp($line);

   if ( ($line =~ /PT in Progress/) ){
      $pt = 1;
      $num_trans++;
      next;
   }

   my ($step, $index) = (split ' ', $line)[0,2]; 

   #if the last set was a PT transfer write output
   if ($pt){
      my $label;
      if ($index == $prev_index){ 
         $label = "Rattle";
         $num_rat++;
      }
      else{
         $label = "PT";
         $num_pt++;
      }

      printf " %9g   -> %9g | %4d   -> %4d  %-s \n", $prev_step, $step, $prev_index, $index, $label;

   }
   
   #Set up the conditions for the next loop
   ($prev_step, $prev_index) = ($step, $index);
   $pt = 0;
}

printf "\nTotal number of steps:        %7d \n", $prev_step;
printf "Total Number of PT:           %7d \n", $num_pt; 
printf "Total Number of Rattles:      %7d \n", $num_rat; 
printf "Total Number of Transitions:  %7d \n\n", $num_trans; 
