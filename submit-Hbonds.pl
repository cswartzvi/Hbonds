#! /usr/bin/perl -w
#

use strict;
use File::Temp qw/ tempfile tempdir /;


my $input_template  = <<EOF;
&system
   filePos     =  '/home/charles/Desktop/Research/Simulations/MD_water-Ions-64/PBE_novdw_330K/OH/Simulation_Files/h2o-64.pos'
   fileCel     =  '/home/charles/Desktop/Research/Simulations/MD_water-Ions-64/PBE_novdw_330K/OH/Simulation_Files/h2o-64.cel'    
   step        =  90000
   ntsp        =  2
   Ospecies    =  1
   Hspecies    =  2
   nsp(1)      =  64
   nsp(2)      =  127
   rcut        =  1.1655
   Hbond_cut   =  3.4
   Hbond_angle =  145.0
/
EOF

while (<>){

   next (/Total/);
   
   #split line to find the proton receiving water step and index and PT/Rattle
   my ($step, $final_step, $index, $type) = (split ' ', $_)[0,2,6,7];

   if ($type ~~ /PT/){
      
      $input_template =~ s/^(step\s+=\s+)\d+$/$1$step/m;
      print $input_template;

      #number of accept donated bonds
      #my @ouput = `echo \"$input_template\" | /home/charles/Desktop/Research/Ext_Programs/HBonds/HBonds.x $index`;
 
      #for (@output){
      #   if (/accepted bonds/){
      #      $num_accept
   } 
}
   
