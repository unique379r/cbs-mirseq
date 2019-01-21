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

# # #=============================================================================================================#
# # NOTE: User must provide all inputs via corresponding input argument file
# # All arguments are compulsory, please DO NOT LEAVE IT EMPTY
# # All arguments are compulsory, please DO NOT LEAVE IT EMPTY
# # This script is a wrapper for module3a and module3b.
# # Module3a Description: Filtering known annotation and prediction of Novel miR
# # Module3b Description: Target Gene prediction of miRNA (might be DE or Novel miR) and their GO, Network analysis
# # Date of modification: 15/4/2016
# # Author: Kesharwani RK email: bioinforupesh2009.au@gmail.com
# # version: 1.0
# # Copyright (c) 2019 Kesharwani RK
# # ================================================================================================================#


# # Usage: bash ./CBS-miRSeq.module3.sh Input_info/Module3_Input.txt


					########################################################
################### WE REQUEST YOU TO PLEASE DO NOT TOUCH BELOW THIS LINE ####################
					########################################################
clear
print_softInfo ()
{
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^Welcome to CBS-miRSeq Module 3^^^^^^^^^^^^^^^^^^^^^^^^#"
echo -e "#		Analysis    : Novel miRNA discovery and target gene prediction	 #"
echo -e "#		Requested Citation: Kesharwani RK et al.(2019)	 #"
echo -e "#		Author		  : bioinforupesh2009.au@gmail.com				 #"
echo -e "#		Copyright (c) 2019 Kesharwani RK				 #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
}
print_softInfo
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

##input file
Module3_Input="$1"
print_USAGE()
{
echo -e "\n"
echo -e "#Please Provide a File Input !!"
echo -e "\n"
echo -e "USAGE:
bash ./CBS-miRSeq.module3.sh Input_Info/Module3_Input.txt"
echo -e "\n"
}
## check file input
if [[ $# != 1 ]]; then
	print_USAGE
	exit 0;
fi
## check right input
if [ ! -s $Module3_Input ]; then
	echo -e "ERROR : Module3_Input=$Module3_Input does not exist\n";
	print_USAGE
	exit 1;
fi

## Extract parameters from the given file
OUTPUT_DIR=$(cat $Module3_Input | grep -w '^OUTPUT_DIR' | cut -d '=' -f2)
MIR_ANNOTATION=$(cat $Module3_Input | grep -w '^MIR_ANNOTATION' | cut -d '=' -f2)
RNAC_ANNOTATION=$(cat $Module3_Input | grep -w '^RNAC_ANNOTATION' | cut -d '=' -f2)
MIR_MATURE_FASTA=$(cat $Module3_Input | grep -w '^MIR_MATURE_FASTA' | cut -d '=' -f2)
REL_MATURE_FASTA=$(cat $Module3_Input | grep -w '^REL_MATURE_FASTA' | cut -d '=' -f2)
PRE_FASTA=$(cat $Module3_Input | grep -w '^PRE_FASTA' | cut -d '=' -f2)
SPS=$(cat $Module3_Input | grep -w '^SPS' | cut -d '=' -f2)
SPECISE_NAME=$(cat $Module3_Input | grep -w '^SPECISE_NAME' | cut -d '=' -f2)
INPUT_QUERY_FASTA=$(cat $Module3_Input | grep -w '^INPUT_QUERY_FASTA' | cut -d '=' -f2)
TARGET_UTR_FASTA=$(cat $Module3_Input | grep -w '^TARGET_UTR_FASTA' | cut -d '=' -f2)
HYBRIDIZATION_THRESHOLD=$(cat $Module3_Input | grep -w '^HYBRIDIZATION_THRESHOLD' | cut -d '=' -f2)
entrezID=$(cat $Module3_Input | grep -w '^entrezID' | cut -d '=' -f2)
ID_Type=$(cat $Module3_Input | grep -w '^ID_Type' | cut -d '=' -f2)
targetHub=$(cat $Module3_Input | grep -w '^targetHub' | cut -d '=' -f2) 
pathways=$(cat $Module3_Input | grep -w '^pathways' | cut -d '=' -f2) 
Ontology=$(cat $Module3_Input | grep -w '^Ontology' | cut -d '=' -f2) 
internet=$(cat $Module3_Input | grep -w '^internet' | cut -d '=' -f2)
plotType=$(cat $Module3_Input | grep -w '^plotType' | cut -d '=' -f2) 
OrgDb=$(cat $Module3_Input | grep -w '^OrgDb' | cut -d '=' -f2)
pvalueCutoff=$(cat $Module3_Input | grep -w '^pvalueCutoff' | cut -d '=' -f2) 
qvalueCutoff=$(cat $Module3_Input | grep -w '^qvalueCutoff' | cut -d '=' -f2) 



## print input from input file
echo -e "\n"
echo -e "#RESULTS DIR:" $OUTPUT_DIR
echo -e "#GENOMIC FEATURES FROM miRBase:" $MIR_ANNOTATION
echo -e "#GENOMIC FEATURES FROM RNAcentral:" $RNAC_ANNOTATION
echo -e "#miRBase REFERENCE SPECIES MATURE FASTA:" $MIR_MATURE_FASTA
echo -e "#RELATED SPECIES MATURE FASTA:" $REL_MATURE_FASTA
echo -e "#miRBase REFERENCE HAIRPIN FASTA:" $PRE_FASTA
echo -e "#REFERENCE SPECIES:" $SPS
echo -e "#NAME OF REFERENCE SPECIES:" $SPECISE_NAME
echo -e "#QUERY_FASTA (Usually DE miRNAs):" $INPUT_QUERY_FASTA
echo -e "#TARGET UTR FASTA:" $TARGET_UTR_FASTA
echo -e "#HYBRIDIZATION THRESHOLD(miRNA:mRNA):" $HYBRIDIZATION_THRESHOLD
echo -e "#UTR fasta contains Entrez gene IDs(yes/no)?:" $entrezID
echo -e "#if NO then what is the first ID (after chromosome) of genes in UTR fasta(SYMBOL/REFSEQ/ENSEMBL)?:" $ID_Type
echo -e "#miRNA:RNA hub to display(small/big/FullNetwork):" $targetHub
echo -e "#Ontology pathway(reactome/kegg):" $pathways
echo -e "#Gene Ontology analysis(DO(for hsa only)/ANY(BP,MF,CC)]):" $Ontology
echo -e "#Internet is accessible(yes/no)?:" $internet
echo -e "#Plot type(pdf/ps)?:" $plotType
echo -e "#Annotation database of your analysis organism:" $OrgDb
echo -e "#P-value cutoff to select genes for GO/pathways output:" $pvalueCutoff
echo -e "#Q-value cutoff to select genes for GO/pathways output:" $qvalueCutoff


## some warnings before to start
echo -e "\n"
echo -e "\t\t#### WARNING !! Before to Run Program ####"
my_color_text "#Make sure internet connection works properly." cyan
echo -e "#Ignore Ineternet connectivity IN CASE you set NO for -> internet"
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


mkdir -p logs
#clear
echo -e '\0033\0143'
print_softInfo () {
echo ''
echo -e "^^^^^^^^^^^^^^^step to run main CBS-miRSeq Modules 3^^^^^^^^^^^^^^^^^^^^^^^^^"
echo -e "		   Analysis:	  Depend on module2a & module2b"
echo -e "	 Requested Citation: Kesharwani RK,Bono E,Chiesa M et al.(2015)"
echo -e "		     Author:  bioinforupesh2009.au@gmail.com#"
echo -e "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
echo ''
}
print_softInfo
echo -e "\n"
echo -e "\t\t\t\t~ ~ ~ ~ ~ CBS-miRSeq.module 3 v1.0 ~ ~ ~ ~ ~"
echo -e "\t\t\t\tAnalysis date:" `date`
echo -e "\n"

# # # # # ### Running analysis by module3a and 3b.

echo -e "\n"
my_color_text "Running Module 3a..." cyan
echo -e "\n"

if [[ -f "./Scripts/CBS-miRSeq.module3a_v1.0.sh" ]]; then
	echo -e "##########################"
	echo -e "## Running Module 3a #####"
	echo -e "##########################"

	bash ./Scripts/CBS-miRSeq.module3a_v1.0.sh "$OUTPUT_DIR/reads_mapped" $OUTPUT_DIR $MIR_ANNOTATION $RNAC_ANNOTATION \
	$OUTPUT_DIR/Index/$SPS*.fa $MIR_MATURE_FASTA $REL_MATURE_FASTA $PRE_FASTA $SPS $SPECISE_NAME
	echo -e "Analysis has been done for Module 3a\n"
else
	echo -e "please make sure you are in the same directory of CBS-miRSeq
	and have provided all necessary input in Module3_Input.txt!!\n"
	exit 1;
fi

sleep 10

echo -e "\n"
my_color_text "Running Module 3b..." cyan
echo -e "\n"

if [[ -f "./Scripts/CBS-miRSeq.module3b_v1.0.sh" ]]; then
	echo -e "##########################"
	echo -e "## Running Module 3b #####"
	echo -e "##########################"

	bash ./Scripts/CBS-miRSeq.module3b_v1.0.sh $INPUT_QUERY_FASTA $TARGET_UTR_FASTA $SPECISE_NAME $SPS $OrgDb $HYBRIDIZATION_THRESHOLD \
	$targetHub $entrezID $ID_Type $pathways $Ontology $pvalueCutoff $qvalueCutoff $internet $plotType $OUTPUT_DIR
else
	echo -e "please make sure you are in the same directory of CBS-miRSeq
	and have provided all necessary input in Module3_Input.txt!!\n"
	exit 1;
fi

set +a

echo -e "#Analysis of Module3 has finished:" `date`
echo -e "#Please go through your Results dir to check your results:" $OUTPUT_DIR
my_color_text "#Thank you for using CBS-miRSeq pipeline." cyan


exit $?

