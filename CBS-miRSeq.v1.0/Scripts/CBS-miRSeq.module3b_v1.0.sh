#!/bin/bash

((
#===============================CBS-miRSeq Module 3 analysis steps ====================================
# Description: Gene Target prediction of miRNA and their GO, Network analysis
# USAGE: bash ./CBS-miRSeq.module3b_v1.0.sh <query_miR_fasta> <UTR_FASTA> <species 3 letter code> <energy_threshold> <output_dir> <entrezID(YES/NO)> <targetHub(small/big/FullNetwork)> <pathways(reactome/kegg)> <internet(yes/no)> <plotType(pdf/eps)> <Ontology(DO/ANY)> <ID_Type(SYMBOL/REFSEQ/ENSEMBL)[##case sensitive]>
# input_fasta = Fasta of your miRNA. it can be known (DE) miRNA or Novel miRNA
# output_dir = directory where results should be written
# 3 PRIME UTR FASTA = 3 PRIME mRNA UTR Fasta of your species. Should be Downloaded from Ensembl biomart.
# sps = 3 letter code of your species
# energy_threshold = A negative value is required for filtering
# Date: 10 April 2015
# version: 1.0
#========================================================================================================

echo -e '\0033\0143'
print_softInfo ()
{
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^CBS-miRSeq Module 3 analysis^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
echo -e "#				~ ~ ~ CBS-miRSeq Module 3b v1  ~ ~ ~				#"
echo -e "#		Analysis    : known and novel miRNA target prediction				#"
echo -e "#		Result: Prediction of Target gene, Gene enrichment, network and pathway analysis#"
echo -e "#		Requested Citation: Kesharwani RK et al.(2019)			#"
echo -e "#		Author			  : bioinforupesh2009.au@gmail.com					#"
echo -e "#		Copyright (c) 2019 Kesharwani RK				 #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
}


# # #####################################################################################
# # # #Example:
# # input_fasta="/data/.../DE_fasta.fa" ## a query fasta in order to find target gene and further to obtain functional enrichment
# # UTR_FASTA="/data/.../genome_3utr.fa"  ## three prime UTR fasta of reference genome (can be download from biomart; refer to CBS-miRSE manual)
# # organism="human" ## Full NAME OF THE species
# # sps="hsa" ## three letter code of reference species
# # OrgDb="org.Hs.eg.db" ## package of annotation Data base of your analysis organism
# # energy_threshold="-20" ## A negative value is required for filtering (i.e -10 or -20 etc) hybridization
# # targetHub="small" ### "small" (recommended) or "big" or "FullNetwork" ## case sensitive
# # entrezID="no" ## (YES; NO)
# # ID_Type="SYMBOL"##ID should be one of: SYMBOL,REFSEQ,ENSEMBL (## case sensitive)
# # pathways="kegg" ## "reactome" or "kegg" ## case sensitive
# # internet="kegg" ## "yes" or "no"
# # plotType="pdf" ### "pdf" or "eps" for plotting
# # Ontology="any" ## gene ontology ## "DO" or "ANY" ; ANY ==> means except DO (i.e. BP,MF,CC)
# # pvalueCutoff="0.05" ## p value cutoff to select genes for GO output
# # qvalueCutoff="0.1" ## q value cutoff to select genes for GO output
# # output_dir="/data/.../target" ## a directory where program can write the results
# # #####################################################################################

print_softInfo
mkdir -p logs
echo -e "\n"
echo -e "\t\t\t~ ~ ~ ~ ~ CBS-miRSeq.module3b_v1.0~ ~ ~ ~ ~"
echo -e "\t\t\t~ ~ ~ ~ ~Target prediction Analysis~ ~ ~ ~ ~"
echo -e "\t\t\tAnalysis date:" `date`

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

## Inputs: command line arguments
input_fasta=$1 ## a query fasta in order to find target gene and further to obtain functional enrichment
UTR_FASTA=$2  ## three prime UTR fasta of reference genome (can be download from biomart; refer to CBS-miRSE manual)
organism=$3 ## Full NAME OF THE species
sps=$4	## three letter code of reference species
OrgDb=$5 ## package of annotation Data base
energy_threshold=$6 ## A negative value is required for filtering
targetHub=$7 ### "small" (recommended) or "big" or "FullNetwork" ## case sensitive
entrezID=$8 ## (YES;NO)
ID_Type=$9 ##ID should be one of: SYMBOL,REFSEQ,ENSEMBL (## case sensitive)
pathways=${10} ## "reactome" or "kegg" ## case sensitive
Ontology=${11} ## gene ontology ## "DO" or "ANY" ; ANY ==> means except DO (i.e. BP,MF,CC)
pvalueCutoff=${12} ## p value cutoff to select genes for GO output
qvalueCutoff=${13} ## q value cutoff to select genes for GO output
internet=${14} ## "yes" or "no"
plotType=${15} ### "pdf" or "eps" for plotting
output_dir=${16} ## a directory where program can write the results


#### sps input checker
ref_cnt=$(echo -n "$sps" | wc -m)
## counts number of charecter in OrgDb
OrgDb_cnt=$(echo -n "$OrgDb" | wc -m)

print_USAGE()
{
	echo -e "USAGE: 
	bash ./CBS-miRSeq.module3b_v1.0.sh <query_miR_fasta> <UTR_FASTA> <organism(NAME OF THE species)> <species 3 letter code> \
<OrgDb(annotationDB like org.Hs.eg.db)> <energy_threshold> <targetHub(small/big/FullNetwork)> <entrezID?(YES/NO)> <GeneID(SYMBOL/REFSEQ/ENSEMBL)[##case sensitive]> \
<pathways(reactome/kegg)> <Ontology(ANY/DO)> <pvalueCutoff> <qvalueCutoff> <internet accessible?(yes/no)> <plotType(pdf/eps)> <output_dir>\n"
	echo -e "EXAMPLE: 
	bash ./CBS-miRSeq.module3b_v1.0.sh /../DE_fasta.fa /../genome_3utr.fa human hsa org.Hs.eg.db -20 small NO SYMBOL kegg DO 0.05 0.1 yes pdf /../Results"
	echo -e "\n"
}


### checking input arguments
if [[ $# -ne 16 ]]; then
	my_color_text "#please supply all inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit;
elif [[ ! -d "$output_dir" ]]; then
    echo -e "Input directory "$output_dir" is not a valid directory !!"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
    exit 1;
elif [[ $ref_cnt -ne 3 ]]; then
	echo -e "\n#Please provide 3 letter code of your reference species !!\n"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ $OrgDb_cnt -ne 12 ]]; then
	echo -e "\n#Probably provided annotation db (OrgDb) is not a valid package name!!\n"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit 1;
else
	echo -e "#job running...\n"
fi

mkdir -p $output_dir/miR_Target
## remove temp files (if present)
rm -f $output_dir/miR_Target/*clean.fa

# echo -e "#fasta file testing..."
# ## input fasta validation
fas=${input_fasta##*/}
if [[ ${fas: -3} == ".fa" ]] || [[ ${fas: -6} == ".fasta" ]]; then
	# # clean invisible carriage return sign from original fasta
	cat $input_fasta | perl -ane 'print "$F[0]\n";' > $output_dir/miR_Target/input_fa_clean.fa
	##echo -e "#microRNA query:" "$fas"
else
	echo -e "ERROR: input fasta file is not a fasta file\n"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
exit 1;
fi

if [[ $ID_Type == "ENSEMBL" ]]; then
	### UTR validation
	if ! grep -q '\|^ENS' "$UTR_FASTA"; then
		echo -e "\nERROR !! \nGiven Target UTR has not appropriate fasta header;
	please refer to CBS-miRSEq manual in order to download appropriate target UTR from biomart.
	or Please make sure header has only 2 or 3 columns
	i.e. ID should be one of: SYMBOL,REFSEQ,ENSEMBL
	Example:
	>chr|ensembl_id|gene_symbol/entrezID
	AGGGCCCGTCCCCAGCCCGGGCCGTCCATCCTC
	ACCATGCCGACCACAAAGGTGTCTGCGGAAACT\n
	or
	>chr|ID
	AGGGCCCGTCCCCAGCCCGGGCCGTCCATCCTC
	ACCATGCCGACCACAAAGGTGTCTGCGGAAACT\n"
	exit 1;
	fi
fi

## UTR cleaning
utr=${UTR_FASTA##*/}

if [[ ${utr: -3} == ".fa" ]] || [[ ${utr: -6} == ".fasta" ]]; then
	# # clean invisible carriage, Sequence unavailable and blank spaces return sign from biomart fasta (mart_export.txt)
	cat $UTR_FASTA | grep -v "Sequence unavailable" | awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' | sed /^$/d | perl -ane 'print "$F[0]\n";' | sed '/^Sequence/d' > $output_dir/miR_Target/UTR_fa_clean.fa
	##echo -e "#3 PRIME UTR Fasta:" "$utr"
else
	echo -e "ERROR: genome UTR file is not a fasta file\n"
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
exit 1;
fi


## check full name of organism
## convert upper to lower case

organismSps=($(echo $organism | tr '[:upper:]' '[:lower:]'))

## convert  lower to upper case of ID if mistaken
ID_Type=($(echo $ID_Type | tr '[:lower:]' '[:upper:]'))

# # energy_threshold validation
en=$(echo -n "$energy_threshold" | wc -c) ## counts numeric character including -ve sign (and it should be 3)
if [[ "$en" -ne 3 ]]; then
	echo -e "ERROR: provided "$energy_threshold" is not a correct value like -10,
please provide 2 digit integers with negative sign."
	my_color_text "#please supply all and correct inputs.." cyan
	echo -e "\n"
	print_USAGE
exit 1;
else
	echo -e "#=============Your Input============#"
	echo -e "#MiRNA fasta:" $input_fasta
	echo -e "#UTR FASTA:" $UTR_FASTA
	echo -e "#Organism:" $organism
	echo -e "#Three letter code of organism:" $sps
	echo -e "#Energy Threshold:" $energy_threshold
	echo -e "#Your UTR fasta contains EntrezID?:" $entrezID
	echo -e "#Gene ID of your UTR fasta:" $ID_Type
	echo -e "#Netwrok hub of miRNA-mRNA:" $targetHub
	echo -e "#Pathways:" $pathways
	echo -e "#Internet:" $internet
	echo -e "#PlotType:" $plotType
	echo -e "#Annotation database of your analysis organism:" $OrgDb
	echo -e "#Type of Ontology[ANY/DO]:" $Ontology
	echo -e "#P-value cutoff to select genes for GO/pathways output:" $pvalueCutoff
	echo -e "#Q-value cutoff to select genes for GO/pathways output:" $qvalueCutoff
	echo -e "#Output directory:" $output_dir
	echo -e "\n"
fi

# # # # ################################# Analysis begins #############################################

###############################
##### RNAhybrid Prediction ####
###############################

### first out put (pair -> miR:Target)
echo -e "#RNAhybrid prediction begins.., please have  patience !! \n"
RNAhybrid -c -f 2,7 -d 1.9,0.28 -p 0.1 -e $energy_threshold -b 1 -t $output_dir/miR_Target/UTR_fa_clean.fa -q $output_dir/miR_Target/input_fa_clean.fa - | tr ':' ' ' | awk '{print $1,"\t", $3}' \
 | sort -u > $output_dir/miR_Target/RNAhybrid_miR_Target_$sps.tab

## checking if program run successfully
if [[ ! -s "$output_dir/miR_Target/RNAhybrid_miR_Target_$sps.tab" ]]; then
		# echo -e "RNAhybrid prediction done.\n"
	# else
		echo -e "Warning !! Output of RNAhybrid is Empty;
there is something wrong either in query miR or in Target UTR fasta
(please examine manually and run again this module)
OR try with low thermodynamics threshold."
	exit 1;
fi


#### Second out put (pair with Structure-> miR:Target)
RNAhybrid -f 2,7 -d 1.9,0.28 -p 0.1 -e $energy_threshold -b 1 -t $output_dir/miR_Target/UTR_fa_clean.fa -q $output_dir/miR_Target/input_fa_clean.fa \
 > $output_dir/miR_Target/RNAhybrid_miR_Target_str_$sps.txt

## checking if program run successfully
 if [[ ! -s $output_dir/miR_Target/RNAhybrid_miR_Target_str_$sps.txt ]]; then
		# echo -e "RNAhybrid prediction done.\n"
	# else
		echo -e "Warning !! Output of RNAhybrid is Empty;
there is something wrong either in query miR or in Target UTR fasta
(please examine manually and run again this module)
OR try with low thermodynamics threshold."
	exit 1;
fi
 echo -e "#RNAhybrid prediction done.\n"

###############################
##### miRanda Prediction ######
###############################


echo -e "#miRanda prediction begins.., please have  patience !! \n"

miranda $output_dir/miR_Target/input_fa_clean.fa $output_dir/miR_Target/UTR_fa_clean.fa -en $energy_threshold -out $output_dir/miR_Target/miranda_miR_Target_$sps.tmp

## checking if program run successfully
if [[ ! -s $output_dir/miR_Target/miranda_miR_Target_$sps.tmp ]]; then
		# echo -e "miRanda prediction done.\n"
	# else
		echo -e "Warning !! Output of RNAhybrid is Empty;
there is wrong something either in query miR or in Target UTR fasta
(please examine manually and run again this module)
OR try with low thermodynamics threshold."
	exit 1;
fi


## Prepare Structure (Hybrid) output
grep -B15 '^>' $output_dir/miR_Target/miranda_miR_Target_$sps.tmp > $output_dir/miR_Target/miranda_miR_Target_str_$sps.txt

grep -i 'Performing Scan:' $output_dir/miR_Target/miranda_miR_Target_str_$sps.txt \
 | sed 's/ /\t/g' | awk '{print $5, "\t", $3}' | sort -u > $output_dir/miR_Target/miranda_miR_Target_$sps.tab

echo -e "#miRanda prediction done.\n"
 # delete temp files
sleep 5
rm -f $output_dir/miR_Target/miranda_miR_Target_$sps.tmp

###################################
### Intersect miRNA:gene Target ###
###################################

## check the format of tab file
## it could be either "chr|EnsembleID|geneSymbol + miR" or "chr|EnsembleID + miR"
tab_count=$(head -1 $output_dir/miR_Target/RNAhybrid_miR_Target_$sps.tab | sed 's/|/\t/g' | awk '{print NF}')

if [[ $tab_count == 4 ]]; then
	## -> parse when chr|EnsembleID|geneSymbol
	cat $output_dir/miR_Target/RNAhybrid_miR_Target_$sps.tab | tr "|" "\t" | sed 's/ /\t/g'\
	| awk '{print $4,"\t", $3}' | sort -u | sed '/^$/d' | grep '[^[:space:]]' > $output_dir/miR_Target/mirna.RNAhybrid.geneTarget.sort

	cat $output_dir/miR_Target/miranda_miR_Target_$sps.tab | tr "|" "\t" | sed 's/ /\t/g'\
	| awk '{print $4,"\t",$3}' | sort -u | sed '/^$/d' | grep '[^[:space:]]' > $output_dir/miR_Target/mirna.miranda.geneTarget.sort
else
	## -> parse when chr|EnsembleID
	cat $output_dir/miR_Target/miranda_miR_Target_$sps.tab | tr "|" "\t"\
	| awk '{print $3,"\t",$2}' | sort -u | sed '/^$/d' | grep '[^[:space:]]' > $output_dir/miR_Target/mirna.miranda.geneTarget.sort

	cat $output_dir/miR_Target/RNAhybrid_miR_Target_$sps.tab | tr "|" "\t"\
	| awk '{print $3,"\t",$2}' | sort -u | sed '/^$/d' | grep '[^[:space:]]' > $output_dir/miR_Target/mirna.RNAhybrid.geneTarget.sort 
fi


echo -e "Preparing Intersect (Common in both algorithm) Target genes...\n"

# # cat $output_dir/miR_Target/RNAhybrid_miR_Target_$sps.tab | tr "|" "\t" | awk '{print $4,"\t", $3}' \
 # # | sort -u | grep '[^[:space:]]' > $output_dir/miR_Target/mirna.RNAhybrid.geneTarget.sort

# # cat $output_dir/miR_Target/miranda_miR_Target_$sps.tab | tr "|" "\t" | awk '{print $4,"\t", $3}' \
 # # | sort -u | sed /^$/d > $output_dir/miR_Target/mirna.miranda.geneTarget.sort

## this was union uniq pair
 ## For GO and Network analysis by R
###paste --delimiter=\\n --serial $output_dir/miR_Target/mirna.RNAhybrid.geneTarget.sort $output_dir/miR_Target/mirna.miranda.geneTarget.sort \
 ##| sort -u | grep '[^[:space:]]' > $output_dir/miR_Target/union.mirna.geneTarget.txt

#### get intersect pair
sort $output_dir/miR_Target/mirna.RNAhybrid.geneTarget.sort $output_dir/miR_Target/mirna.miranda.geneTarget.sort | uniq -d > $output_dir/miR_Target/intersect.genes.uniq.pair.txt

## check final results and clean temp dir
if [[ -f $output_dir/miR_Target/intersect.genes.uniq.pair.txt ]]; then
	## make header
	sed -i '1d' $output_dir/miR_Target/intersect.genes.uniq.pair.txt
	sed -i -e '1imiRNA\tmRNA' $output_dir/miR_Target/intersect.genes.uniq.pair.txt
	echo -e "#Selection of Overlapped (Intersect) genes, done.\n"
	## remove temp files
	rm -f $output_dir/miR_Target/*.sort 2>/dev/null
	rm -f $output_dir/miR_Target/*clean.fa 2>/dev/null
	echo -e "\n_miRNA Target Analysis finished": `date`
else
	sleep 5
	## remove temp files
	rm -f $output_dir/miR_Target/*.sort 2>/dev/null
	rm -f $output_dir/miR_Target/*clean.fa 2>/dev/null
fi


##NOTE: echo -e "Please use generated "intersect.genes.uniq.pair.txt"
##files to analyse GO and Network analysis respectively by R Module 3 of CBS-miRSeq pipeline\n"

############################ calling R Script for GO and Network analysis #######################

### Important Note: might one error occur due to bugs in Packages, but NO worries, results will be fine.
### Probable error: [1] "ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method."

echo -e "#Running GO enrichment and Network analysis..\n"

## Calling Rscript for GO analysis..
GO_script="./Scripts/CBS-miRSeq.GO_analysis.R"

if [[ -s $output_dir/miR_Target/intersect.genes.uniq.pair.txt ]]; then

	if [[ -f $GO_script ]]; then

		Rscript $GO_script $output_dir/miR_Target/intersect.genes.uniq.pair.txt $output_dir/miR_Target $organism $sps $OrgDb $entrezID $targetHub $pathways $internet $plotType $Ontology $ID_Type $pvalueCutoff $qvalueCutoff

	else
		echo -e "CBS-miRSeq.GO_analysis.R script not found !! please make sure that you are in
CBS-miRSeq script directory\n"
	fi
else
	echo -e "mirna.geneTarget lists are empty !!
GO and Network analysis halted.\n"
	exit;
fi


# # echo -e Important Note: might one error occur due to bugs in R Packages, but NO worries, results will be fine.
# # echo -e something like: "#ERROR: The estimated pi0 <= 0......,please ignore this error"

echo -e "\n######"
echo -e "_Gene Target prediction of miRNA and their GO, Network Analysis finished\n": `date`
echo -e "######"

) 2>&1) | tee -a module3b.log

if [[ -e $logs/module3b.log ]]; then

	rm -f $logs/module3b.log

fi



mv -f module3b.log logs

echo -e "\n#log file also has created into the directory of CBS-miRSeq: logs/module3b.log\n"


exit $?

