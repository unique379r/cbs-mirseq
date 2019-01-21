#!/bin/bash -xv

((
#=======================================CBS-miRSeq Module 2b ==================================================
# Description: isomiR detection from mapped
# USAGE: ./CBS-miRSeq.module2b_v1.0.sh <mapped/bam_dir> <ref.genomeFASTA> <mirna_annotation/gff3> <known_mirnaFASTA> <sps 3 letter code> <output_dir># mapped_reads/bam = mapping reads directory, aligned by CBS-miRSeq.module 1
# mapped/bam_dir = mapping directory contained sorted bam files
# ref.genomeFASTA = directory of genome FASTA
# mirna_annotation = miRBase gff3 file of your sps
# mature_known_mirnaFASTA=mature fasta of your sps (can be downloaded from miRBase ftp) ## file
# sps_3_letter_code= three lketter code of your species ## string
# output_dir = directory where results should be written
# Note: Script should be run in script folder of CBS-miRSeq
# Date : 14 April 2015
# Authors: Keshrwani RK email: bioinforupesh2009.au@gmail.com
#==========================================================================================================

################################################################
# # Example:
# # mapped/bam_dir="/data/..../reads_mapped" ## directory
# # ref.genomeFASTA="/data/...../Index/genome.fa" ## file
# # mirna_annotation/gff3="/data/.../hsa.gff3" ## file
# # mature_known_mirnaFASTA="/data/.../hsa.mature.fa" ## file
# # sps_3_letter_code="hsa" ## string
# # output_dir="/data/..../Results" ## directory
###############################################################

# clear
echo -e '\0033\0143'
print_softInfo ()
{
echo -e "#^^^^^^^^^^^^^^^^^^^^CBS-miRSeq Module 2 detection of isomiR^^^^^^^^^^^^^^^^^^^^#"
echo -e "#				~ ~ ~ CBS-miRSeq Module 2 v1.0	~ ~ ~		#"
echo -e "#		Analysis    : iso-miR detection					#"
echo -e "#		Result: iso-miR html results corresponds to known mirna		#"
echo -e "#		Requested Citation: Kesharwani RK et al.(2019)#"
echo -e "#		Author			  : bioinforupesh2009.au@gmail.com			#"
echo -e "#		Copyright (c) 2019 Kesharwani RK				 #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
}

print_softInfo

echo -e "\n"
echo -e "\t\t\t~ ~ ~ ~ ~ CBS-miRSeq.module2b_v1.0~ ~ ~ ~ ~"
echo -e "\t\t\t~ ~ ~ ~ ~iso-miR detection~ ~ ~ ~ ~"
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
bam="$1"
ref_genome="$2"
mirna_annotation="$3"
mirna_fa="$4"
sps="$5"
out="$6"

print_USAGE()
{
	echo -e "USAGE: 
	bash ./CBS-miRSeq.module2b_v1.0.sh <reads_mapped_dir> <ref.genomeFASTA> <mirna_annotation/gff3> <All.sps.mature.faFASTA> <sps 3 letter code> <output_dir>\n"
	echo "EXAMPLE: 
	bash ./CBS-miRSeq.module2b_v1.0.sh /../reads_mapped /../Index/genome.fa /../annotation/hsa.gff3 /../annotation/All.sps.mature.fa hsa /../Results"
	echo -e "\n"
}

ref_cnt=$(echo -n "$sps" | wc -m)

if [[ $# -ne 6 ]]; then
	my_color_text "#please supply all inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit;
elif [[ $ref_cnt -ne 3 ]]; then
	echo -e "\nPlease provide 3 letter code of your reference species !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ ! -d "$bam" ]] && [[ ! -d "out" ]]; then
echo -e "\nalignment directory: "$bam"
#OR
Output directory: "$out"
are not a valid directory, please check !!"
my_color_text "#please supply all and correct inputs" cyan
echo -e "\n"
print_USAGE
else
	echo -e "#Testing files for isomiR detection proceeded..\n"
fi

echo -e "#species :" "$sps"
### testing input files
# genome validation
gen=${ref_genome##*/}
if [[ ${gen: -3} == ".fa" ]] || [[ ${gen: -6} == ".fasta" ]]; then
	echo -e "#Genome:" "$gen"
else
	echo -e "ERROR: genome file is not a fasta file\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
print_USAGE
	exit 1;
fi
# # mature miR validation
miR=${mirna_fa##*/}
if [[ ${miR: -3} == ".fa" ]] || [[ ${miR: -6} == ".fasta" ]]; then
	echo -e "#miR fasta:" "$miR"
else
	echo -e "ERROR: miRBase file is not a fasta file\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
print_USAGE
	exit 1;
fi

# # miRBase annotation file validation
miR_gff3=${mirna_annotation##*/}
if [[ ${miR_gff3: -5} == ".gff3" ]]; then
	echo -e "#miR annotation:" "$miR_gff3"
else
	echo -e "ERROR: miRBase annotation is not a valid gff3 file format !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
print_USAGE
	exit 1;
fi

echo -e "#annotation file testing done, analysis proceeded.\n"
 # # gff3 to gff2 and extract fasta as hairpin by adding 35 flanking
 # # extract hairpin from gff3 (gff3 to gff2)

if [[ -f $mirna_annotation ]]
then
	##echo -e "gff3 found..\n"
	## echo $mirna_annotation ## original files and path by user
	miR_dir=$(dirname "${mirna_annotation}") ## extract path only
	## echo "$miR_dir"
	## base-name
	##echo "${mirna_annotation##*/}"
	touch $miR_dir/$sps.gff3.gff2
	awk 'NR==4;NR==5;NR==6' $mirna_annotation >> $miR_dir/$sps.gff3.gff2
	echo -e "# File type: gff2 hairpin" >> $miR_dir/$sps.gff3.gff2
	echo -e "# `date`" >> $miR_dir/$sps.gff3.gff2
	echo -e "#" >> $miR_dir/$sps.gff3.gff2
	### its globally work for all sps (indeed its not needed to be sort for other sps)
	grep -v 'Zv9\|^#' $mirna_annotation | sed -n '/miRNA_primary_transcript/p' | sort -k1n	\
	| sed 's/Alias=MI[[:digit:]]\{7\};//g;s/=/="/g;s/;/";/g;s/$/";/g;s/miRNA_primary_transcript/miRNA/g' \
	>> $miR_dir/$sps.gff3.gff2
	echo -e "gff3 to gff2 conversion done !!\n"
else
	echo -e "miRBase annotation .gff3 does not found, please provide a valid directory of input files !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
print_USAGE
	exit 1;
fi

# #  make gff2 hairpin to 35flank+fasta from genome
# #  samtool must be install

if [[ -f $ref_genome ]] && [[ -f ./Scripts/ExtractGenomeSequences.pl ]]
then
	echo -e "sequence extraction for gff2 being processed...please wait\n"
	## v1.1
	perl ./Scripts/ExtractGenomeSequences.pl -gff $miR_dir/$sps.gff3.gff2 -fa $ref_genome -s $sps -flank 35 > $miR_dir/$sps.35nt.pre.fa
	echo -e "done\n"
else
	echo -e "genome fasta does not found or you are not inside the CBS-miRSeq scripts directory !! \n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
print_USAGE
	exit 1;
fi

mkdir -p $out/isomiRs

# # ## miRSpring call
# # # #	 Note: bam must be sorted and have their bai index
map=($bam/*.sorted.bam)
#bai=($bam/*.sorted.bam.bai)
if [[ -f "${map[0]}" ]] && [[ -f $mirna_fa ]]
then
	echo -e "bam found..\n"
	for i in "${map[@]}"
	do
		echo -e "$ isomiR detection proceeded..\n"
		## v1.3
		perl ./Scripts/BAM_to_Intermediate.pl -ml 0 -s $sps -flank 35 -pre $miR_dir/$sps.35nt.pre.fa -gff $miR_dir/$sps.gff3.gff2 -mat $mirna_fa -bam $i -ref 0 -out $out/isomiRs/${i##*/}.miRspring.input.txt
		## v1.2
		perl ./Scripts/Intermediate_to_miRspring.pl -flank 35 -in $out/isomiRs/${i##*/}.miRspring.input.txt -s $sps -out $out/isomiRs/${i##*/}.miRspring.html
		echo -e "isomiR analysis done !!
check your results..\n"
	done
else
	echo -e "bam or mature fasta is not found..!!\n"
	exit 1;
fi

rm -f $miR_dir/*.gff2
rm -f $miR_dir/*35nt.pre.fa

echo -e "\nAnalysis finished": `date`

) 2>&1) | tee -a module2b.log

if [[ -s $logs/module2b.log ]]; then

	rm -f $logs/module2b.log

fi

mv -f module2b.log logs
echo -e "\n#log file also has created into the directory of CBS-miRSeq: logs/module2b.log\n"

exit $?

