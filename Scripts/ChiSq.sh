#!/bin/bash

#script to calculate the chi squared function for the fitted parameters

#outline first
#take in parameter range, and experimental parameters, feed into fdjacobian (outputs theoretical cross sections at x_i)
#first create grid, will need 3 loops
#calculate total chi squared
#figure out the best way to formate this file with these numbers

#echo "sigma_1: "
#read sig1
#echo "sigma_2: "
#read sig2
#echo "sigma_3: "
#read sig3

#sig1=$(echo "scale=15; $sig1*2" | bc)
#sig2=$(echo "scale=15; $sig2*2" | bc)
#sig3=$(echo "scale=15; $sig3*2" | bc)

sig1=2954.424280552161
sig2=12.9303345442257
sig3=3.2493168069258633
sta1=$(echo "scale=15; 131.18593407820873-2*$sig1" | bc)
sta2=$(echo "scale=15; 0.93402076452647564-2*$sig2" | bc)
sta3=$(echo "scale=15; 0.96835868174987560-2*$sig3" | bc)
end1=$(echo "scale=15; 131.18593407820873+2*$sig1" | bc)
end2=$(echo "scale=15; 0.93402076452647564+2*$sig2" | bc)
end3=$(echo "scale=15; 0.96835868174987560+2*$sig3" | bc)

#while [ $cou1 -lt $end1 ]; do
#   while [ $cou2 -lt $end2 ]; do
      while [$cou3 -lt  ] ; do
         echo "$sta1 $sta2 $sta3" >  parm.txt
#	 /mnt/analysis/fresco-errors/bin/fdjacobian -fd_h 0.001 -x_in parm.txt
	 sta3=$(echo "scale=15; $sta3+1" | bc)
	 echo "$sta3"
      done 

#eventually have it read this in from the Jacobian
