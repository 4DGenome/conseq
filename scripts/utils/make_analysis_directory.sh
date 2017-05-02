#!/bin/bash


#==================================================================================================
# Created on: 2017-04-06
# Usage: ./make_analysis_directory.sh
# Author: javier.quilez@crg.eu
# Goal: makes directory for analysis within the project directory
#==================================================================================================

CONSEQ=/users/GR/mb/jquilez/projects/conseq

# variables
project=$1
analysis=$2

# check variables are passed as script parameters
if [ -n "$project" ] && [ -n "$analysis" ]; then
	my_date=`date +"%Y-%m-%d"`
	directory_name=${my_date}_${analysis}
else
	echo -e "\nusage: make_analysis_directory.sh <project> <analysis>\n"
	exit
fi

# make directories and files
ANALYSIS=$CONSEQ/projects/$project/analysis/${my_date}_${analysis}
mkdir -p $ANALYSIS/{scripts,tables,figures,data}
md=$ANALYSIS/$directory_name.md
rm -f $md
echo "# $directory_name" >> $md
echo `printf '%100s\n' | tr ' ' -` >> $md
echo -e "\n**objective: ...**" >> $md
echo -e "\n**paths are relative to $CONSEQ**\n\n" >> $md

echo -e "\nanalysis directory created at $ANALYSIS\n"
