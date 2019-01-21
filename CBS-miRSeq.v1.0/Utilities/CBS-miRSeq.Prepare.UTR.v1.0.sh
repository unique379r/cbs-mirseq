#!/bin/bash


## Date of creation: 21/07/2016
## This script is a part of the CBS-miRSeq help to prepare UTR which must downloaded from biomart
## Author:  Kesharwani RK; bioinforupesh2009.au@gmail.com

echo -e '\0033\0143'
print_softInfo () {
echo ''
echo -e "#^^^^^^^^^^^^^^^Welcome to the Utility script for UTR prepration^^^^^^^^^^^^^#"
echo -e "#		   Aim:	  To Prepare UTR		     	     	     #"
echo -e "#	 Requested Citation: Kesharwani RK et al.(2019)	     #"
echo -e "#		     Author:  Kesharwani RK; bioinforupesh2009.au@gmail.com	     #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
echo ''
}
print_softInfo

echo -e "\t\t\tAnalysis date:" `date`


#Inputs: command line arguments
Output_dir="$1" ## OutPut directory
sps="$2" ## three letter code of reference species
biomartfasta="$3" ## input

## for checking length of sps
ref_cnt=$(echo -n "$sps" | wc -m)

print_USAGE()
{
echo -e "USAGE: bash ./CBS-miRSeq.Prepare.UTR.v1.0.sh <Output_dir> <reference_species_3_letter_code> <biomartfasta/mart_export.txt>\n"
}

# checking input arguments
if [[ $# -ne 3 ]]; then
	echo -e "#Please supply all inputs.."
	echo -e "\n"
	print_USAGE
	echo -e "#====================================#"
	echo -e "NOTE: Probable reference species:"
	echo -e "mmu <= Mouse"
	echo -e "hsa <= Human"
	echo -e "dre <= Zebrafish"
	echo -e "rno <= Rat"
	echo -e "#=====================================#"
	echo -e "\n"
	exit 0;
elif [[ ! -d "$Output_dir" ]]; then
    echo -e "\n#Output directory "$Output_dir" \ndoes not exist !!\n"
	echo -e "#please supply all and correct inputs"
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ $ref_cnt -ne 3 ]]; then
	echo -e "\n#Please provide 3 letter code of your reference species !!\n"
	echo -e "\nError(2)... \nplease supply all inputs.. \n"
	print_USAGE
	exit 1;
fi


## Analysis
## Count given fasta including "Sequence unavailable"
echo -e "#No. of given fasta (Included Sequence unavailable)..\n"
grep -c '>' $biomartfasta
echo -e "\n"
echo -e "Processing....\n"
## Remove Sequence unavailable
cat $biomartfasta | grep -v "Sequence unavailable" | awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' | sed /^$/d | perl -ane 'print "$F[0]\n";' | sed '/^Sequence/d' > $Output_dir/sps.UTR_fa_clean.fa
## Clean invisible carriage
perl -ane 'print "$F[0]\n";' $Output_dir/sps.UTR_fa_clean.fa > $Output_dir/$sps.UTR_fa_clean.fa
## paste ENS id where no symbol to maintain the required format >1|ENSDARG00987654321|Genesymbol
sed -i -e 's/|ENSDARG[[:digit:]]\{11\}$/&&/g' $Output_dir/$sps.UTR_fa_clean.fa

## remove temp dir
rm -rf $Output_dir/sps.UTR_fa_clean.fa

## Count fasta with sequences
echo -e "#No. of Actual fasta..\n"
grep -c '>' $Output_dir/$sps.UTR_fa_clean.fa
echo -e "Output name of the file:" $Output_dir/$sps.UTR_fa_clean.fa
echo -e "\n"
echo -e "#Your UTR fasta has been Prepared, please check and launch your analysis using module 3 of CBS-miRSeq.\n"
echo -e "#Thank you for using the CBS-miRSeq.\n"

