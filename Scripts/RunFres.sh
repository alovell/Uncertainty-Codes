#!/bin/bash

#script to run fresco multiple times based on changing set of fitted parameters

var=( 0 0 0 0 0 0 0 0 0 0 0 0 0 )
#lines to find in fresco input (will have to be changed depending on case)
#rv1
#in1="&POT kp=1 type=1 p1=131.13 p2=0.93421 p3=0.96832  p4=1.7 p5=1.32 p6=0.35  /"
#in2="&POT kp=1 type=2 p1=0 p2=0 p3=0  p4=1.7 p5=1.32 p6=0.35  /"
#in3="&POT kp=1 type=3 p1=10.6 p2=1.36 p3=0.85   /"
#wv1
#in1="&POT kp=1 type=1 p1=95.910000 p2=1.150000 p3=0.780000 p4=1.700000 p5=1.320000 p6=0.350000  /"
#in2="&POT kp=1 type=2 p1=0 p2=0 p3=0 p4=1.700000 p5=1.320000 p6=0.350000  /"
#in3="&POT kp=1 type=3 p1=10.600000 p2=1.360000 p3=0.850000   /"
#wv2
in1="&POT kp=1 type=1 p1=124.68 p2=1.1318 p3=0.8853  p4=0.0 p5=1.35 p6=0.73  /"
in2="&POT kp=1 type=2 p1=0 p2=0 p3=0  p4=13.1445 p5=1.200 p6=0.549775  /"
in3="&POT kp=1 type=3 p1=7.8 p2=1.16 p3=0.84   /"
#fresco input file
iFile="dwba-fit2"
file="dwba-fit2.in"
#parameters file
input=("$@")
fFile=${input[0]}

#replace variables with numbers
sed -i "s/aso /12/g" $fFile
sed -i "s/rso /11/g" $fFile
sed -i "s/Vso /10/g" $fFile
sed -i "s/Wsa /9/g" $fFile
sed -i "s/Wsr /8/g" $fFile
sed -i "s/Ws /7/g" $fFile
sed -i "s/Wwa /6/g" $fFile
sed -i "s/Wwr /5/g" $fFile
sed -i "s/Wv /4/g" $fFile
sed -i "s/a /3/g" $fFile
sed -i "s/r /2/g" $fFile
sed -i "s/V / 1/g" $fFile



while read line
do
   if [[ $line == "#"* ]];
   then
     #figure out which parameters were being fit
     #put non-zeros in var array for fitted parameters
     IFS=' ' read -a parm <<< "$line"
     #echo "${parm[1]}"
     for i in "${!parm[@]}"
     do 
       #echo "$i ${parm[i]}"
       if [[ $i != 0 ]]; then
       n=${parm[$i]}
#       echo "$n"
       var[$n]=$n
       fi
     done
#     for j in "${!var[@]}"
#     do 
#       echo "$j ${var[$j]}"
#     done

   else 
     #read in parameters to the correct places and then run fresco
     IFS=' ' read -a v0 <<< "$line"
     #need seds and correct lines for them
     #parse lines and put into var array if value is zero.
     IFS=' ' read -a v1 <<< "$in1"
     IFS=' ' read -a v2 <<< "$in2"
     IFS=' ' read -a v3 <<< "$in3"
     m=0
     #puts values that change into var
     for x in {1..12}
     do
       if [[ ${var[$x]} != 0 ]]; then 
#         echo "var start ${var[$x]}"
         var[$x]=${v0[$m]}
#	 echo "m $m var after ${var[$x]}"
	 m=$m+1
#       echo "$x $line"
#	 want=${line#*\ }
#	 echo "$want"
#	 var[$x]="$want"
#	 line=${line%\ *}
#	 echo "$line"
       fi
     done
     
     #puts values from 1st line that don't change into var
     for x in {3..8}
     do
       if [[ ${var[$x-2]} == 0 ]]; then
         p=${v1[$x]}
	 #echo "$p"
         var[$x-2]=${p#*=}
       fi
     done
     
     #puts values from 2nd line that don't change into var
     for x in {6..8}
     do
       if [[ ${var[$x+1]} == 0 ]]; then
         q=${v2[$x]}
	 #echo "$q"
         var[$x+1]=${q#*=}
       fi
     done
     
     #puts values from 3rd line that don't change into var
     for x in {3..5}
     do
       if [[ ${var[$x+7]} == 0 ]]; then
         r=${v3[$x]}
	 #echo "$r"
         var[$x+7]=${r#*=}
       fi
     done
     
     #writes out var array to make sure everything's in the right place
     #COMMENT OUT when script is working
     for y in "${!var[@]}"
     do 
       echo "$y ${var[$y]}"
     done
     
     #replacing variables using awk
     count=$(cat ${file} | wc -l)
     #echo "There are $count lines"
     
     #create temp input file
     temp="${iFile}_TEMP"
     rm -f $temp
     
     for ((n=1;n<=$count;n++)) do
     current=$(head -$n $file | tail -1)
     testbit=$(echo $current | awk '{if ($0 ~/&POT kp=1 type=1/){print 1}else if ($0 ~/&POT kp=1 type=2/){print 2}else if ($0 ~/&POT kp=1 type=3/){print 3}else{print 0}}')
     
     if [[ $testbit -eq 0 ]]; then
       echo "$current" >> $temp
     elif [ $testbit -eq 1 ]; then 
#       echo "$testbit"
       printf " &POT kp=1 type=1 p1=%f p2=%f p3=%f p4=%f p5=%f p6=%f /\n" ${var[1]} ${var[2]} ${var[3]} ${var[4]} ${var[5]} ${var[6]} >> $temp
     elif [ $testbit -eq 2 ]; then
       printf " &POT kp=1 type=2 p1=0 p2=0 p3=0 p4=%f p5=%f p6=%f /\n" ${var[7]} ${var[8]} ${var[9]} >> $temp
     else 
       printf " &POT kp=1 type=3 p1=%f p2=%f p3=%f /\n" ${var[10]} ${var[11]} ${var[12]} >> $temp
     fi
     done
     
     cp -f $temp ./$file
     rm -f $temp

     #run fresco
     ~/fresco/fresco < $file > $iFile.out
     
     #then probably need to copy all output/input/fort into a new directory - talk to Filomena about how much should be saved
     
   fi 
   
   #very bulk-ily change any var spot back to zero if that variable was not fitted with sfresco
   for z in {0..12}
   do
     for w in {1..12}
     do 
       if [[ $z != ${parm[1]} ]] && [[ $z != ${parm[2]} ]] && [[ $z != ${parm[3]} ]] && [[ $z != ${parm[4]} ]] && [[ $z != ${parm[5]} ]] && [[ $z != ${parm[6]} ]] && [[ $z != ${parm[7]} ]] && [[ $z != ${parm[8]} ]] && [[ $z != ${parm[9]} ]] && [[ $z != ${parm[10]} ]] && [[ $z != ${parm[11]} ]] && [[ $z != ${parm[12]} ]]; then
#         echo "$z ${parm[w]} yay!"
         var[$z]=0
#	 echo "$z ${var[$z]}"
       else
         var[$z]=${var[$z]}
#         echo "$z ${var[$z]}"
       fi
     done
   done 
done < $fFile


#list of possible parameters (usual 12 that we would fit)
#V,r,a, Wv,Wwr,Wwa
#Ws,Wsr,Wsa
#Vso,rso,aso

#this script would have to be modified if the real surface term were fitted
