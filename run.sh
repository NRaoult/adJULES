#!/usr/bin/env bash
set -eu #e=stop at first error, u=fail undefinted variable, x=print lines 

here=$(pwd)

#if [ "$3" = "choose" ]
#then  
#echo "Input the set of sites you wish to use (seperated by commas):"
#read choice
#set -- "${@:1:2}" "$choice" "${@:4}"
#echo "$@"
#fi 

#if [ "$5" = "on" ]
#then  
#echo "Pick nudge factor:"
#read choice
#set -- "${@:1:5}" "$choice"
#echo "$@"
#fi 

Rscript create_temp.R "$@"

dir="temp/temp_$1"
cd $dir
Rscript $here/$dir/run.R > $here/$1_output.out 2>&1 &
 
rm -fr $dir
