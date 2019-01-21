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

# # =========================================================================================#
# # NOTE: User must provide all inputs via corresponding input argument file
# # All arguments are compulsory, please DO NOT LEAVE IT EMPTY
# # This script is a wrapper for module2a and module2b.
# # Module2a Description: Differential analysis and classification of Ensembl biotypes
# # Module2b Description: isomiR prediction in every sample
# # Date of modification: 15/4/2016
# # Author: Kesharwani RK email: bioinforupesh2009.au@gmail.com
# # version: 1.0
# # Copyright (c) 2019 Kesharwani RK
# # ==========================================================================================#
# # Note:
# # Please make sure that your miRBase matrix (hsa.miRBase.report.txt) is in same order such as
# # Control Control Control Treatment Treatment Treatment....
# # Healthy Healthy Healthy Unhealthy Unhealthy Unhealthy....
# # WildType WildType WildType GeneKnockOut GeneKnockOut GeneKnockOut....
# # In general----> groupA vs groupB
# # This script is a wrapper for module2a and module2b.
# # Module2a Description: Differential Expression analysis
# # Module2b Description: isomiR detection from mapped sam
# # =========================================================================================#

# # Usage: bash ./CBS-miRSeq.module2.sh Input_info/Module2_Input.txt


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
echo -e "#^^^^^^^^^^^^^^^Welcome to CBS-miRSeq Module 2^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
echo -e "#		Analysis  : DE, Ensembl biotype classifications & iso-miR detection#"
echo -e "#		Requested Citation: Kesharwani RK et al.(2019)	   #"
echo -e "#		Author		  : bioinforupesh2009.au@gmail.com				   #"
echo -e "#		Copyright (c) 2019 Kesharwani RK				   #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
}
print_softInfo

##input file
Module2_Input="$1"
print_USAGE()
{
echo -e "\n"
echo -e "#Please Provide a File Input !!"
echo -e "\n"
echo -e "USAGE:
bash ./CBS-miRSeq.module2.sh Input_Info/Module2_Input.txt"
echo -e "\n"
}
## check file input
if [[ $# != 1 ]]; then
	print_USAGE
	exit 0;
fi
## check right input
if [ ! -s $Module2_Input ]; then
	echo -e "ERROR : Module2_Input=$Module2_Input does not exist\n";
	print_USAGE
	exit 1;
fi

## Extract parameters from the given file
OUTPUT_DIR=$(cat $Module2_Input | grep -w '^OUTPUT_DIR' | cut -d '=' -f2)
MIR_ANNOTATION=$(cat $Module2_Input | grep -w '^MIR_ANNOTATION' | cut -d '=' -f2)
ALL_SPS_MATURE_FASTA=$(cat $Module2_Input | grep -w '^ALL_SPS_MATURE_FASTA' | cut -d '=' -f2)
MIR_MATURE_FASTA=$(cat $Module2_Input | grep -w '^MIR_MATURE_FASTA' | cut -d '=' -f2)
SPS=$(cat $Module2_Input | grep -w '^SPS' | cut -d '=' -f2)
conditionA=$(cat $Module2_Input | grep -w '^conditionA' | cut -d '=' -f2)
startCol_grpA=$(cat $Module2_Input | grep -w '^startCol_grpA' | cut -d '=' -f2)
conditionB=$(cat $Module2_Input | grep -w '^conditionB' | cut -d '=' -f2)
startCol_grpB=$(cat $Module2_Input | grep -w '^startCol_grpB' | cut -d '=' -f2)
LowCountsfeaturesFilterByCPM=$(cat $Module2_Input | grep -w '^LowCountsfeaturesFilterByCPM' | cut -d '=' -f2)
CutoffLowCountsfeaturesFilter=$(cat $Module2_Input | grep -w '^CutoffLowCountsfeaturesFilter' | cut -d '=' -f2)
RemoveInconsistentFeatures=$(cat $Module2_Input | grep -w '^RemoveInconsistentFeatures' | cut -d '=' -f2)
ThresholdToRemoveInconsistentFeatures=$(cat $Module2_Input | grep -w '^ThresholdToRemoveInconsistentFeatures' | cut -d '=' -f2)
performIndependentFilter=$(cat $Module2_Input | grep -w '^performIndependentFilter' | cut -d '=' -f2)
plotType=$(cat $Module2_Input | grep -w '^plotType' | cut -d '=' -f2)
pval_Cutoff=$(cat $Module2_Input | grep -w '^pval_Cutoff' | cut -d '=' -f2)
padj_Cutoff=$(cat $Module2_Input | grep -w '^padj_Cutoff' | cut -d '=' -f2)
log2FC_Cutoff=$(cat $Module2_Input | grep -w '^log2FC_Cutoff' | cut -d '=' -f2)


## print input from input file
echo -e "\n"
echo -e "#RESULTS DIR:" $OUTPUT_DIR
echo -e "#GENOMIC FEATURES FROM miRBase:" $MIR_ANNOTATION
echo -e "#miRBase REF SPECIES MATURE FASTA:" $MIR_MATURE_FASTA
echo -e "#miRBase ALL SPECIES MATURE FASTA:" $ALL_SPS_MATURE_FASTA
echo -e "#REFERENCE SPECIES:" $SPS
echo -e "#condition/GroupA:" $conditionA
echo -e "#GroupA.startCol in miRBase Expression Matrix:" $startCol_grpA
echo -e "#condtiion/GroupB:" $conditionB
echo -e "#GroupB.startCol in miRBase Expression Matrix:" $startCol_grpB
echo -e "#Filter low/unexpressed counts using CPM(yes/no[In case of filtration by raw counts])?:" $LowCountsfeaturesFilterByCPM
echo -e "#CutoffLowCountsfeaturesFilter:" $CutoffLowCountsfeaturesFilter
echo -e "#Remove hypervariants from Matrix(yes/no)?:" $RemoveInconsistentFeatures
echo -e "#Threshold (CV cutoff) to remove hypervariants:" $ThresholdToRemoveInconsistentFeatures
echo -e "#Perform Independent Filtering?(yes/no)" $performIndependentFilter
echo -e "#Plot type(pdf/ps)?:" $plotType
echo -e "#P-Value Cutoff:" $pval_Cutoff
echo -e "#P-Adj(FDR) Cutoff:" $padj_Cutoff
echo -e "#Log2Fold Change Cutoff:" $log2FC_Cutoff


## some warnings before to start
echo -e "\n"
echo -e "\t\t#### WARNING !! Before to Run Program ####"
my_color_text "#Make sure internet connection works properly." cyan
echo -e "#Ignore Ineternet connectivity IN CASE, you already have installed all required packages using Rscript ./CBS-miRSeq.Required.Packages.R\n"
my_color_text "#Please Recheck Your Arguments Printed Above." red
echo -n "Continue ? (y/n) : "
read ans
if [[ "${ans}" != "y" ]] && [[ "${ans}" != "Y" ]]; then
	clear
	echo -e "\n"
	my_color_text "Please note that any missing inputs/arguments can cause to fail your analysis" red
	my_color_text "Please recheck arguments inside Module2_Input.txt and Run this script again." cyan
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
echo -e "^^^^^^^^^^^^^^^step to run main CBS-miRSeq Modules 2^^^^^^^^^^^^^^^^^^^^^^^^^"
echo -e "		   Analysis:	  Depend on module2a & module2b"
echo -e "	 Requested Citation: Kesharwani RK,Bono E,Chiesa M et al.(2015)"
echo -e "		     Author:  bioinforupesh2009.au@gmail.com"
echo -e "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
echo ''
}
print_softInfo
echo -e "\n"
echo -e "\t\t~ ~ ~ ~ ~ CBS-miRSeq.module 2 v1.0 ~ ~ ~ ~ ~"
echo -e "\t\tAnalysis date:" `date`
echo -e "\n"

## check input matrix

if [[ -s $OUTPUT_DIR/quantification/$SPS.ensembl.report.txt && -s $OUTPUT_DIR/quantification/$SPS.miRBase.report.txt ]]; then
	echo -e "\n"
	echo -e "#Quantified Ensembl Biotype:" $OUTPUT_DIR/quantification/$SPS.ensembl.report.txt
	echo -e "#Quantified miRBase miRNAs:" $OUTPUT_DIR/quantification/$SPS.miRBase.report.txt
	echo -e "\n"
else
	my_color_text "ERROR!! Expression matrices are empty" red
	echo -e "#Please check quantification directory in:" $OUTPUT_DIR
	exit 1;
fi

### Running analysis by module2a and 2b.

echo -e "Running Module 2a..\n"

if [[ -f "./Scripts/CBS-miRSeq.module2a_v1.0.sh" ]]; then
	echo -e "##########################"
	echo -e "## Running Module 2a #####"
	echo -e "##########################"

	bash ./Scripts/CBS-miRSeq.module2a_v1.0.sh $OUTPUT_DIR/quantification/$SPS.ensembl.report.txt $OUTPUT_DIR/quantification/$SPS.miRBase.report.txt \
	$OUTPUT_DIR $SPS $conditionA $startCol_grpA $conditionB $startCol_grpB $MIR_MATURE_FASTA \
	$LowCountsfeaturesFilterByCPM $CutoffLowCountsfeaturesFilter $RemoveInconsistentFeatures $ThresholdToRemoveInconsistentFeatures $performIndependentFilter $plotType \
	$pval_Cutoff $padj_Cutoff $log2FC_Cutoff
	echo -e "Analysis of Module 2a finished.\n"
else
	echo -e "please make sure you are in same the directory of CBS-miRSeq
	and have provided all necessary input in Module2_Input.txt!!\n"
	exit 1;
fi

sleep 10

echo -e "Running Module 2b..\n"

 if [[ -f "./Scripts/CBS-miRSeq.module2b_v1.0.sh" ]]; then
	 echo -e "##########################"
	 echo -e "## Running Module 2b #####"
	echo -e "##########################"

	 bash ./Scripts/CBS-miRSeq.module2b_v1.0.sh $OUTPUT_DIR/reads_mapped $OUTPUT_DIR/Index/$SPS*.fa $MIR_ANNOTATION $ALL_SPS_MATURE_FASTA $SPS $OUTPUT_DIR
	 echo -e "Analysis of Module 2b finished\n"
 else
	 echo -e "please make sure you are in same the directory of CBS-miRSeq
	 and have provided all necessary input in Module2_Input.txt!!\n"
	 exit 1;
 fi

set +a

echo -e "#Analysis of Module2 has finished:" `date`
echo -e "#Now you may proceed with Module 3\n"
echo -e "#Please go through your Results dir to check your results:" $OUTPUT_DIR
my_color_text "#Thank you for using CBS-miRSeq pipeline." cyan

exit $?

