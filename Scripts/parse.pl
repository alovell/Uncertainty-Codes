#!/usr/bin/perl
use strict;
use warnings;

#parsing chi_file to make contour plots

my $iFile = $ARGV[0];
my $oFile = "";
my $num1;
my $num2;
my @bestfit;

#open chi_file
open (FILE, '<', $iFile) or die "Could not open $iFile:  $!";

#test
my $i=1;
#read in file and separate chunks between "headers"
while (<FILE>) {
    #get rid of new line
    chomp;
    #split line based on white spaces, put into line array
    my @line = split(/ +/);
    if ($i==0) {
      while (my ($it,$el) = each @line) {
        $bestfit[$it]=$el;
        #print "$bestfit[$it] \n";
      }
    }
    else {
      if ($line[-1]==2 || $line[-1]==3 || $line[-1]==4 || $line[-1]==5 || $line[-1]==6 || $line[-1]==7 || $line[-1]==8 || $line[-1]==9 || $line[-1]==10 || $line[-1]==11 || $line[-1]==12) {
        $num1=$line[1];
	$num2=$line[-1];
	#print "$num1$num2 \n";
	$oFile = 'fit3_90Zrpgrid' . $num1 . $num2 . '.txt';
	print "$oFile \n";
	open(TEMP, '>>', $oFile) or die "Could not open output file:  $!";
      }
      else {
        print TEMP "$line[$num1]  $line[$num2]  $line[4] \n";
	#print TEMP "$line[$num1]  $line[$num2] \n";
      }
    }
    $i=1
    
}