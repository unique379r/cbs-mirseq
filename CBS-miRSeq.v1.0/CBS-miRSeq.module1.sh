#!/bin/bash

set -e
set -a

# # This file is part of CBS-miRSeq.

# # CBS-miRSeq is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.

# # CBS-miRSeq is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.

# # You should have received a copy of the GNU General Public License
# # along with CBS-miRSeq.  If not, see <http://www.gnu.org/licenses/>

# # ====================================================================================================================#
# # NOTE: User must provide all inputs via corresponding input argument file
# # All arguments are compulsory, please DO NOT LEAVE IT EMPTY
# # This script is a wrapper for module1a and module1b.
# # Module1a Description:  building.index of species genome as well download genome if not present (require internet connectivity if genome is not already downloaded)
# # Module1b Description: QC (1st and 2nd),Trimming,Mapping,Quantification and summarization of datasets
# # Date of modification: 15/4/2016
# # Author: Kesharwani RK email: bioinforupesh2009.au@gmail.com
# # version: 1.0
# # Copyright (c) 2019 Kesharwani RK
# # ====================================================================================================================#

# # Usage: bash ./CBS-miRSeq.module1.sh Input_info/Module1_Input.txt

					########################################################
################### WE REQUEST YOU TO PLEASE DO NOT TOUCH BELOW THIS LINE ####################
					########################################################
## function for color
my_color_text(){
    local exp=$1;
    local color=$2;
    if ! [[ $color =~ '^[0-9]$' ]]; then
		case $(echo $color | tr '[:upper:]' '[:lower:]') in
		black) color=0 ;;
		red) color=1 ;;
		green) color=2 ;;
		yellow) color=3 ;;
		blue) color=4 ;;
		magenta) color=5 ;;
		cyan) color=6 ;;
		white|*) color=8 ;; ## white or invalid color
		esac
    fi
    tput setaf $color;
    echo $exp;
    tput sgr0;
}
clear
print_softInfo ()
{
echo -e "#^^^^^^^^^^^^^^^^^^^^Welcome to CBS-miRSeq Module 1^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
echo -e "#		Analysis    : building.ref.index, QC, Trimming, Mapping and Quantification#"
echo -e "#		Requested Citation: Kesharwani RK et al.(2019)		  #"
echo -e "#		Author		  : bioinforupesh2009.au@gmail.com					  #"
echo -e "#		Copyright (c) 2019 Kesharwani RK				 	  #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
}
print_softInfo
##input file
Module1_Input="$1"
print_USAGE()
{
echo -e "\n"
echo -e "#Please Provide a File Input !!"
echo -e "\n"
echo -e "USAGE:
bash ./CBS-miRSeq.module1.sh Input_info/Module1_Input.txt"
echo -e "\n"
}
## check file input
if [[ $# != 1 ]]; then
	print_USAGE
	exit 0;
fi
## check right input
if [ ! -s $Module1_Input ]; then
	echo -e "ERROR : Module1_Input=$module1_Input does not exist\n";
	print_USAGE
	exit 1;
fi

## Extract parameters from the given file

DIR_REF_GENOME=$(cat $Module1_Input | grep -w '^DIR_REF_GENOME' | cut -d '=' -f2)
SPS=$(cat $Module1_Input | grep -w '^SPS' | cut -d '=' -f2)
INDEX_COLOR=$(cat $Module1_Input | grep -w '^INDEX_COLOR' | cut -d '=' -f2)
ENS_RELEASE=$(cat $Module1_Input | grep -w '^ENS_RELEASE' | cut -d '=' -f2)
ENS_AssemblyGRCh=$(cat $Module1_Input | grep -w '^ENS_AssemblyGRCh' | cut -d '=' -f2)
DIR_INPUT_READS=$(cat $Module1_Input | grep -w '^DIR_INPUT_READS' | cut -d '=' -f2)
MIR_ANNOTATION=$(cat $Module1_Input | grep -w '^MIR_ANNOTATION' | cut -d '=' -f2)
ENS_ANNOTATION=$(cat $Module1_Input | grep -w '^ENS_ANNOTATION' | cut -d '=' -f2)
MY_ADAPTER=$(cat $Module1_Input | grep -w '^MY_ADAPTER' | cut -d '=' -f2)
OUTPUT_DIR=$(cat $Module1_Input | grep -w '^OUTPUT_DIR' | cut -d '=' -f2)
NumberOfSamples=$(cat $Module1_Input | grep -w '^NumberOfSamples' | cut -d '=' -f2)
## print input from input file
echo -e "\n"
echo -e "#DIR OF REFERENCE GENOME:" $DIR_REF_GENOME
echo -e "#REFERENCE SPECIES:" $SPS
echo -e "#BUILD INDEX IN COLOR?:" $INDEX_COLOR
echo -e "#ENSEMBL RELEASE:" $ENS_RELEASE
echo -e "#GENOME ASSEMBLY:" $ENS_AssemblyGRCh
echo -e "#DIR OF INPUT READS:" $DIR_INPUT_READS
echo -e "#NUMBER OF SAMPLE:" $NumberOfSamples
echo -e "#GENOMIC FEATURES FROM miRBase:" $MIR_ANNOTATION
echo -e "#GENOMIC FEATURES FROM ENSEMBL:" $ENS_ANNOTATION
echo -e "#ADAPTER SEQUENCE:" $MY_ADAPTER
echo -e "#RESULTS DIR:" $OUTPUT_DIR

## some warnings before to start
echo -e "\n"
echo -e "\t\t#### WARNING !! Before to Run Program ####"
my_color_text "#Make sure internet connection works properly." cyan
echo -e "#Ignore Ineternet connectivity IN CASE you already have downloaded your reference genome.\n"
my_color_text "#Please Recheck Your Arguments Printed Above." red
echo -n "Continue ? (y/n) : "
read ans
if [[ "${ans}" != "y" ]] && [[ "${ans}" != "Y" ]]; then
	clear
	echo -e "\n"
	my_color_text "Please note that any missing inputs/arguments can cause to fail your analysis" red
	my_color_text "Please recheck arguments inside Module1_Input.txt and Run this script again." cyan
	my_color_text "Thank you for using CBS-miRSeq pipeline." blue
	echo -e "\n"
	exit 0;
fi

########################
#### Analysis begins ###
########################
mkdir -p logs
#clear
echo -e '\0033\0143'
print_softInfo () {
echo ''
echo -e "^^^^^^^^^^^^^^^step to run main CBS-miRSeq Modules 1^^^^^^^^^^^^^^^^^^^^^^^^^"
echo -e "		   Analysis:	  Depend on module1a & module1b"
echo -e "	 Requested Citation: Kesharwani RK,Bono E,Chiesa M et al.(2015)"
echo -e "		     Author:  bioinforupesh2009.au@gmail.com"
echo -e "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
echo ''
}
print_softInfo
echo -e "\n"
echo -e "\t\t  ~ ~ ~ ~ CBS-miRSeq.module 1 v1.0 ~ ~ ~ ~	 "
echo -e "\t\tAnalysis date:" `date`
echo -e "\n"

# # # Running analysis by module1a and 1b.

echo -e "Running Module 1a...\n"

if [[ -f "./Scripts/CBS-miRSeq.module1a_v1.0.sh" ]]; then
	echo -e "##########################"
	echo -e "## Running Module 1a #####"
	echo -e "##########################"

	bash ./Scripts/CBS-miRSeq.module1a_v1.0.sh $DIR_REF_GENOME $SPS $INDEX_COLOR $ENS_RELEASE $ENS_AssemblyGRCh
	echo -e "#Analysis of Module 1a finished.\n"
else
	echo -e "please make sure you are in the same directory of CBS-miRSeq
	and have provided all necessary input in Module1_Input.txt!!\n"
	exit 1;
fi

sleep 10

echo -e "Running Module 1b...\n"

## check if genome.fa and index are present

fasta=$(ls $DIR_REF_GENOME/Index/*.fa 2> /dev/null)

if [[ ! -s $fasta ]]; then
	my_color_text "ERROR!! Genome fasta is empty, 
	please check your input of version and Assembly of Ensembl database." red
	echo -e "Visit->: ftp://ftp.ensembl.org/pub/ \n"
	my_color_text "#Please clean your directory (Except Annotation one; if you are in same dir) 
	and launch this module again." cyan
	exit 1;
else
	## copying dir of Index from Genome to Output dir
	cp -Rf $DIR_REF_GENOME/Index $OUTPUT_DIR
fi

if [[ -f "./Scripts/CBS-miRSeq.module1b_v1.0.sh" ]]; then
	echo -e "##########################"
	echo -e "## Running Module 1b #####"
	echo -e "##########################"

	bash ./Scripts/CBS-miRSeq.module1b_v1.0.sh $DIR_INPUT_READS $OUTPUT_DIR $OUTPUT_DIR/Index $MIR_ANNOTATION $ENS_ANNOTATION $SPS $INDEX_COLOR $ENS_RELEASE $MY_ADAPTER $NumberOfSamples
	echo -e "Analysis of Module 1b finished.\n"
else
	echo -e "please make sure you are in same directory of CBS-miRSeq
	and have provided all necessary input in Module1_Input.txt!!\n"
	exit 1;
fi

set +a


echo -e "#Analysis of Module1 has Finished:" `date`
echo -e "#Now you may proceed with Module 2\n"
echo -e "#Please go through your Results dir to check your results:" $OUTPUT_DIR
my_color_text "#Thank you for using CBS-miRSeq pipeline." cyan

exit $?

