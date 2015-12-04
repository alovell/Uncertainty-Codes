#!/bin/bash

#script to read through all of the input files in the folder with a certain name

m=1
for i in X0*.txt 
do 
   sh /mnt/analysis/fresco-errors/Scripts/RunFres.sh $i
   file1="fort$m.201"
   file2="fort$m.202"
   cp fort.201 ./$file1
   cp fort.202 ./$file2
   let m=m+1
   echo "$file1 $file2"
done

#do sh /mnt/analysis/fresco-errors/Scripts/RunFres $i

#try to read in from command line as a test
#args=("$@")
#echo ${args[0]}
