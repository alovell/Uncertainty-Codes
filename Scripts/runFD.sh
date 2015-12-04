#!/bin/bash

#essentially a wrapper for chi^2 calculation
#run fdjacobian on a grid to get all of the chisquare values 

#counter - 0 first time to calculate grid and sig values
i=0
j=1

if [ $i==0 ]
then 
   #run fdjacobian, calc_chi
   ../bin/fdjacobian -fd_h 0.001 -search_file search-102.in -x_in x.in > out.out
   chical
   #make grid of points 
   rm -r parm_grid.txt
   grid
   #calculat sigma values
   i=1
fi 

#put two lines from grid into x.in
#unless it's a header line
while read line
do
   IFS=' ' read -a array <<< "$line"
   if [ ${array[1]} == 2 ] || [ ${array[1]} == 3 ] || [ ${array[1]} == 4 ] || \
          [ ${array[1]} == 5 ] || [ ${array[1]} == 6 ]; then 
      #IFS=' ' read -a array <<< "$line"
      echo "$line" >> chicompphys.txt
      #echo ${array[1]}
   else
      #if [ $j == 1 ]
      #then
         #echo $j
         echo -n "$line" > x.in
	# j=2
      #else
      #   echo -n "$line " >> x.in
         #echo $j
       #  j=1
      #fi
      #compute jacobian for each set of parameters
      ../bin/fdjacobian -fd_h 0.001 -search_file search-102.in -x_in x.in > out.out
      #calculate chisq value and put it into chi_file with parameters
      chical
   fi
   #echo $line
done < parm_grid.txt

#run fdjacobian
#fdjacobian -fd_h 0.001 -search_file search-102.in -x_in x.in
#calculate the chisq value
#chical 