#!/bin/bash -xv

((
######################## CBS-miRSeq Module 2a ############################
#	Description:Diff Expr and ensemble gene classification script wrapper
#	option should be equal to Rscript
#	Parameters:
#	ensembl_counts= A count matrix of Ensembl biotypes (*.txt)
#	raw_exprs= miRNA Expression matrix obtained by ht-seq/featuresCounts (*.txt)
#	output_dir_path= Path where results should be produced
#	sps= Analysis sps 3 letter code
#	group.A = experimental condition A
#	startCol_grpA= Group First start column number in your miRNA count matrix (usually 2)
#	group.B = experimental condition B
#	startCol_grpB= Group Second start column number in your miRNA count matrix
#   mature_fasta= A fasta file of your analysis species
#	LowCountsfeaturesFilterByCPM= Filter low/unexpressed tags before analysis (yes or no; [default Nothing])
#	CutoffLowCountsfeaturesFilter= [default 1] ->filteration CutoffLowCountsfeaturesFilter either used for LowCountsfeaturesFilterByCPM or Raw counts filter if LowCountsfeaturesFilterByCPM is set to no
#	RemoveInconsistentFeatures= Remove outliers from the matrix before to perform DE [yes or no; default Nothing/Optional]
#	ThresholdToRemoveInconsistentFeatures= [default 1.5] ->Less than this ThresholdToRemoveInconsistentFeatures tags will be keep(applied in both model)
#	performIndependentFilter=[default no] ->Independent filtering if yes then applied to both model
#   plotType= pdf or ps (postscript type)
# 	pval_Cutoff=pvalue cutoff for plotting ## Numeric Value
# 	padj_Cutoff=FDR cutoff for plotting ## Numeric Value
# 	log2FC_Cutoff=lof2fc cutoff for plotting ## Numeric Value
#	date: 15 April 2016
#	version : v1.0
#	Authors: Keshrwani RK email: bioinforupesh2009.au@gmail.com
##############################################################################

###===========================================================================
# # Example:
# # ensembl_counts="/..../quantification/hsa.ensembl.report.txt" ## file
# # raw_exprs="/...../quantification/hsa.miRBase.report.txt" ## file
# # output_dir_path="/..../Results" ## absolute path/dir
# # sps="hsa" ## string
# # conditionA="Control"
# # startCol_grpA="2" ## Numeric value
# # conditionB="Treatment"
# # startCol_grpB="5" ## Numeric Value
# # mature_fasta="/...../hsa.mature.fa" ## file
# # LowCountsfeaturesFilterByCPM="yes" ## string
# # CutoffLowCountsfeaturesFilter="1" ## Numeric value
# # RemoveInconsistentFeatures="yes" ## string
# # ThresholdToRemoveInconsistentFeatures="1.5"
# # performIndependentFilter="no" ## string
# # plotType="pdf" ## or ps[default pdf]
# # pval_Cutoff=pvalue cutoff for plotting ## Numeric Value
# # padj_Cutoff=FDR cutoff for plotting ## Numeric Value
# # log2FC_Cutoff=lof2fc cutoff for plotting ## Numeric Value
###===========================================================================

#clear
echo -e '\0033\0143'
print_softInfo ()
{
echo -e "#^^^^^CBS-miRSeq Module 2 DE and ensembl gene classification^^^^^^^^^^^^^^^^^^#"
echo -e "#				~ ~ ~ CBS-miRSeq Module 2 v1.0	~ ~ ~	      #"
echo -e "#		Analysis    : DE and ensembl biotype classification	      #"
echo -e "#		Result:		DE excel file  and plots		      #"
echo -e "#		Requested Citation: Kesharwani RK et al.(2019)#"
echo -e "#		Author			  : bioinforupesh2009.au@gmail.com		      #"
echo -e "#		Copyright (c) 2019 Kesharwani RK				 #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
}


print_softInfo
echo -e "\n"
echo -e "\t\t\t~ ~ ~ ~ ~ CBS-miRSeq.module2a_v1.0~ ~ ~ ~ ~"
echo -e "\t\t\t~ ~ ~ ~ ~DE and ensembl gene classification~ ~ ~ ~ ~"
echo -e "\t\t\tAnalysis date:" `date`
mkdir -p logs

## function for color
my_color_text(){
    local exp=$1;
    local color=$2;
    if ! [[ $color =~ '^[0-9]$' ]] ; then
       case $(echo $color | tr '[:upper:]' '[:lower:]') in
	black) color=0 ;;
	red) color=1 ;;
	green) color=2 ;;
	yellow) color=3 ;;
	blue) color=4 ;;
	magenta) color=5 ;;
	cyan) color=6 ;;
	white|*) color=7 ;; ## white or invalid color
       esac
    fi
    tput setaf $color;
    echo $exp;
    tput sgr0;
}


## Inputs...
ensembl_counts="$1"
raw_exprs="$2"
output_dir_path="$3"
sps="$4"
conditionA="$5"
startCol_grpA="$6"
conditionB="$7"
startCol_grpB="$8"
mature_fasta="$9"
LowCountsfeaturesFilterByCPM="${10}"
CutoffLowCountsfeaturesFilter="${11}"
RemoveInconsistentFeatures="${12}"
ThresholdToRemoveInconsistentFeatures="${13}"
performIndependentFilter="${14}"
plotType="${15}"
pval_Cutoff="${16}"
padj_Cutoff="${17}"
log2FC_Cutoff="${18}"

## sps input checker
#ref_cnt=$(echo "$sps" | wc -L)
ref_cnt=$(echo -n "$sps" | wc -m)

print_USAGE()
{
	echo -e "USAGE:
	bash ./CBS-miRSeq.module2a_v1.0.sh <ensembl_counts.txt> <miR_expression.txt> <output_dir_path> <analysis sps 3 letter code> \
<condition/group.A> <startCol.grpsA(from miRBase table)> \
<condition/group.B> <startCol.grpsB(from miRBase table)> \
<mature fasta of analysis sps> \
<LowCountsfeaturesFilterByCPM(yes/no)> <CutoffLowCountsfeaturesFilter> \
<RemoveOutliersFromMatrix(yes/no)> <ThresholdToRemoveInconsistentFeatures> \
<performIndependentFilter(yes/no)> <plotType(pdf/eps)> <pval_Cutoff> <pval_Cutoff> <log2FC_Cutoff>\n"
	echo -e "EXAMPLE:
	bash ./CBS-miRSeq.module2a_v1.0.sh /..../quantification/hsa.ensembl.report.txt /...../quantification/hsa.miRBase.report.txt /..../Results hsa Control 2 Treatment 5 /...../hsa.mature.fa no yes 1 yes 1.5 no pdf 0.05 0.05 1"
	echo -e "\n"
}

##checking input arguments
if  [[ $# -ne 18 ]]; then
	echo -e "\nAnalysis date:" `date`
	my_color_text "#please supply all inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit;
elif [[ ! -d "$output_dir_path" ]]; then
	echo -e "\nError... output directory "$output_dir_path" not a valid directory !!"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
    exit 1;
elif [ "$ref_cnt" -ne 3 ]; then
	echo -e "\nError... Please provide 3 letter code of your reference species !!\n"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit 1;
else
	echo -e "job is running...\n"
	echo species : "$sps"
fi

# # mature miR validation
miR=${mature_fasta##*/}
if [[ ${miR: -3} == ".fa" ]] || [[ ${miR: -6} == ".fasta" ]]; then
	echo -e "#miR fasta:" "$miR"
else
	echo -e "ERROR: mature fasta file is a not valid fasta \n"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit 1;
fi

## expression file validation
ens=${ensembl_counts##*/}
if [[ ${ens: -4} == ".txt" ]]; then
	echo -e "#ensembl Expression:" "$ens"
else
	echo -e "ERROR: Ensembl file is a not valid expression file \n"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit 1;
fi

mirEx=${raw_exprs##*/}
if [[ ${mirEx: -4} == ".txt" ]]; then
	echo -e "#miR Expression:" "$mirEx\n"
else
	echo -e "ERROR: miR count file is a not valid expression file \n"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit 1;
fi

mkdir -p $output_dir_path/DE
## Calling Rscript for diff expr analysis..
DE_script="./Scripts/CBS-miRSeq.DE.HTSFilter.script_v1.R"

if [[ -f $DE_script ]]; then

	Rscript $DE_script $raw_exprs $output_dir_path/DE $sps $conditionA $startCol_grpA $conditionB $startCol_grpB \
	$LowCountsfeaturesFilterByCPM $CutoffLowCountsfeaturesFilter $RemoveInconsistentFeatures $ThresholdToRemoveInconsistentFeatures \
	$performIndependentFilter $plotType $pval_Cutoff $padj_Cutoff $log2FC_Cutoff
else
	#echo -e "\nInputs are ok but..."
	echo -e "Error(2)..."
	echo -e "BS-miRSeq.DE.HTSFilter.script_v1.R does not exist in the same directory."
fi

## Calling Rscript for ensembl gene annotation counts table..
Ens_script="./Scripts/CBS-miRSeq.ensembl.gene_v1.R"

if [[ -f $Ens_script ]]; then

	Rscript $Ens_script $ensembl_counts $output_dir_path/DE
else
	#echo -e "\nInputs are ok but..."
	echo -e "Error(3)..."
	echo -e "CBS-miRSeq.ensembl.gene_v1.R does not exist in the same directory."
fi


#### extract Intersect miRNs fasta and convert RNA into DNA

if ls $output_dir_path/DE/*Intersect.merged.stat*.csv 1> /dev/null 2>&1; then

	column -s, -t $output_dir_path/DE/*Intersect.merged.stat*.csv | awk '{print $1}' | sed '1d' > $output_dir_path/DE/DE_miR.txt

	grep -w -A1 -f $output_dir_path/DE/DE_miR.txt $mature_fasta | sed 's/--//g' | sed /^$/d | sed '/^[^>]/ y/uU/tT/' | sed '/^[^>]/ s/n//g' > $output_dir_path/DE/DE_mature.fa ## can be used for target prediction

else
	echo -e "Intersect.merged.stat is empty !!"
	rm -f $output_dir_path/DE/DE_miR.txt
fi

echo -e "\nAnalysis finished": `date`

) 2>&1) | tee -a module2a.log

if [[ -s $logs/module2a.log ]]; then

	rm -f $logs/module2a.log

fi

mv -f module2a.log logs

echo -e "\n#log file also has created into the directory of CBS-miRSeq: logs/module2a.log\n"

exit $?

