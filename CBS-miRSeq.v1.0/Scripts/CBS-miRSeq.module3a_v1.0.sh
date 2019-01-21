#!/bin/bash

((
#=============================CBS-miRSeq Module 3a analysis steps ======================================================
# Description: Filtering known annotation and Novel miR prediction from aligned_reads
# USAGE: ./CBS-miRSeq.module3_v1.sh <input_map> <output_path> <miRBase annotation/GFF3> <RNAcentral/bed> <genomeFASTA> <mature_miR_fasta> <related_sps_mature_fasta> <your_sps_precursor_fasta> <species 3 letter code> <species Full Name>
# input_map= Destination of sam/bam
# output_dir = Directory where results should be written
# miR_anno/GTF = Ensembl annotation in gff (genome co-ordinates) http://www.ensembl.org/info/data/ftp/index.html
# RNAcentral/bed= RNAcentral annotation in bed format (i.e. Homo_sapiens.GRCh38.bed) ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/5.0/genome_coordinates/
# genome= A fasta File where ref genome is located (usually in Index folder (.../Index/hsa.genome_v79.fa)
# mir_fasta = mature fasta file of your sps
# related_sps_fasta= related known mature miRNA fasta of your species
# pre_fasta= miR precursor seq fasta of your sps
# sps = Analysis species 3 letter code (i.e. hsa/dre/mmu/rat..)
# species=Name of species; i.e. hsa=Human or dre=Zebrafish or mmu=Mouse.. #Hyperlink to the UCSC browser entry
# Date : 14 April 2015
# version: 1.0
# Authors: Keshrwani RK email: bioinforupesh2009.au@gmail.com
#==================================================================================================================

#clear
echo -e '\0033\0143'
print_softInfo ()
{
echo -e "#^^^^^^^^^^^^^^^CBS-miRSeq Module 3 Novel miRNA prediction^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
echo -e "#				~ ~ ~ CBS-miRSeq Module 3a v1.0 ~ ~ ~			 #"
echo -e "#		Analysis    : Novel miRNA prediction					 #"
echo -e "#		Result: Novel miRNA prediction by miRDeep2				 #"
echo -e "#		Requested Citation: Kesharwani RK et al.(2019)		 #"
echo -e "#		Author			  : bioinforupesh2009.au@gmail.com				 #"
echo -e "#		Copyright (c) 2019 Kesharwani RK				 #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
}

# # ########################################################################################################
# # #Example: used without quotes in terminal
# # input_map="/data/.../reads_mapped"	## destination of aligned.sam/bam
# # output_dir="/data/.../"
# # miR_anno="/data/.../miRBase/dre.gff3"
# # RNAcentral="/data/.../RNAcentral/Danio_rerio.Zv9.bed"
# # genome="/data/.../dre.genome.fa" ## genome of your species
# # mir_fasta="/data/.../miRBase_v21_mature.fa" ## known mature miRNA fasta of your species
# # related_sps_fasta="/data/.../tni_mature_v21.fa" ## related known mature miRNA fasta of your species
# # pre_fasta="/data/.../dre_pre_v21.fa" ## known precursor miRNA fasta of your species
# # sps="dre" #Analysis species 3 letter code (i.e. hsa/dre/mmu/rat..)
# # species="Zebrafish" #Hyperlink to the UCSC browser entry
# # # # ####################################################################################################

print_softInfo

echo -e "\n"
echo -e "\t\t\t~ ~ ~ ~ ~ CBS-miRSeq.module3a_v1.0~ ~ ~ ~ ~"
echo -e "\t~ ~ ~ ~ ~Novel miRNA prediction and Counts matrix & Novel bed generation~ ~ ~ ~ ~"
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

### Inputs: command line arguments
input_map=$1
output_dir=$2
miR_anno=$3
RNAcentral=$4
genome=$5
mir_fasta=$6
related_sps_fasta=$7
pre_fasta=$8
sps=$9
species=${10}  ## used for miRDeep2 #Hyperlink to the UCSC browser entry


## sps input checker
ref_cnt=$(echo -n "$sps" | wc -c)

print_USAGE()
{
	echo -e "USAGE: 
	bash ./CBS-miRSeq.module3_v1.sh <reads_mapped_dir> <output_path> <miRBase annotation/GFF3> <RNAcentral/bed> <genomeFASTA> <mature_miR_fasta> <related_sps_mature_fasta> <your_sps_precursor_fasta> <species 3 letter code> <Full Name of the species>\n"
	echo -e "EXAMPLE: 
	bash ./CBS-miRSeq.module3a_v1.0.sh /../reads_mapped /../Results /../annotation/hsa.gff3 /../annotation/Homo_sapiens.GRCh38.RNAcentral.bed /../Index/hsa.genome_v84.fa /../annotation/hsa.mature_v21.fa /../annotation/mmu.Related.mature_v21.fa /../annotation/hsa.precursor_v21.fa hsa human"
	echo -e "\n"
}

# # checking input arguments
if [[ $# -ne 10 ]]; then
	my_color_text "#please supply all inputs.." cyan
	echo -e "\n"
	print_USAGE
	exit;
elif [[ ! -d "$input_map" ]] && [[ ! -d "$output_dir" ]]; then
    echo -e "\nError...\n
Input directory: "$input_map"
#OR
Output directory: "$output_dir"
are not a valid directory, please Recheck !!"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ ! -f "$miR_anno" ]] && [[ ! -f "$RNAcentral" ]]; then
	echo -e "\nError...
Either "$miR_anno" OR "$RNAcentral"
is not a valid annotation or missing as input
#Please provide valid annotations  !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ ! -f "$mir_fasta" ]] && [[ ! -f "$related_sps_fasta" ]] && [[ ! -f "$pre_fasta" ]] && [[ ! -f "$genome" ]]; then
	echo -e "\nError...
Either "$genome" OR "$mir_fasta" OR "$related_sps_fasta" OR "$pre_fasta"
is not a valid fasta or missing as input
#Please provide valid fasta  !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ $ref_cnt -ne 3 ]]; then
	echo -e "\nError...
#Please provide 3 letter code of your reference species !!"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
else
	echo -e "\njob running...\n"
fi

echo -e "File Testing...\n"
echo -e "#species :" "$sps"
echo -e "#species Full Name :" "$species"

## declare an array variable
declare -a bam=($input_map/*.bam)

if [[ ! -f "${bam[0]}" ]]; then
	echo -e "mapped bam files not found..
please provide valid directory..\n"
	exit 1;
		else
	rm -rf $input_map/mirdeep_input 2> /dev/null
	mkdir -p $input_map/mirdeep_input
	echo -e "#mapping bam found"
fi

# # checking fasta provided by user
# # clean invisible carriage return sign from original fasta

gen=${genome##*/}
miR=${mir_fasta##*/}
related_miR=${related_sps_fasta##*/}
pre_mir=${pre_fasta##*/}

# genome validation
if [[ ! -f "$input_map/mirdeep_input/genome_clean.fa" ]]; then
	if [[ ${gen: -3} == ".fa" ]] || [[ ${gen: -6} == ".fasta" ]]; then
		cat $genome | perl -ane 'print "$F[0]\n";' > $input_map/mirdeep_input/genome_clean.fa
		echo -e "#Genome fasta:" "$gen"
	else
		echo -e "ERROR: genome file is not a fasta file\n"
		my_color_text "#please supply all and correct inputs" cyan
		echo -e "\n"
		print_USAGE
		exit 1;
	fi
fi
# # mature miR validation
if [[ ! -f "$input_map/mirdeep_input/miR_clean.fa" ]]; then
	if [[ ${miR: -3} == ".fa" ]] || [[ ${miR: -6} == ".fasta" ]]; then
		cat $mir_fasta | perl -ane 'print "$F[0]\n";' > $input_map/mirdeep_input/miR_clean.fa
		echo -e "#mature miRNA fasta:" "$miR"
	else
		echo -e "ERROR: species mature miR file is not a fasta file\n"
		my_color_text "#please supply all and correct inputs" cyan
		echo -e "\n"
		print_USAGE
		exit 1;
	fi
fi

# # mature related sps miR validation
if [[ ! -f "$input_map/mirdeep_input/related_miR_clean.fa" ]]; then
	if [[ ${related_miR: -3} == ".fa" ]] || [[ ${related_miR: -6} == ".fasta" ]]; then
		cat $related_sps_fasta | perl -ane 'print "$F[0]\n";' > $input_map/mirdeep_input/related_miR_clean.fa
		echo -e "#related sps mature miRNA fasta:" "$related_miR"
	else
		echo -e "ERROR: species related_miR file is not a fasta file\n"
		my_color_text "#please supply all and correct inputs" cyan
		echo -e "\n"
		print_USAGE
		exit 1;
	fi
fi

# # pre_mir validation
if [[ ! -f "$input_map/mirdeep_input/pre_mir_clean.fa" ]]; then
	if [[ ${pre_mir: -3} == ".fa" ]] || [[ ${pre_mir: -6} == ".fasta" ]]; then
		cat $pre_fasta | perl -ane 'print "$F[0]\n";' > $input_map/mirdeep_input/pre_mir_clean.fa
		echo -e "#precursor miRNA fasta:" "$pre_mir"
	else
		echo -e "ERROR: species related_miR file is not a fasta file\n"
		my_color_text "#please supply all and correct inputs" cyan
		echo -e "\n"
		print_USAGE
		exit 1;
	fi
else
	echo -e "#fasta file testing done, analysis proceeded\n"
fi

# # miRBase and RNAcentral annotations file validation

miR_gff3=${miR_anno##*/}
RNA_bed=${RNAcentral##*/}
# miR_gff3 validation
if [[ ${miR_gff3: -5} != ".gff3" ]]; then
	echo -e "ERROR: miRBase annotation is not a valid gff3 file format !!\n"
echo -e "\nplease supply all & correct inputs !!
\nUSAGE: ./CBS-miRSeq.module3_v1.sh <input_map> <output_path> <miRBase annotation/GFF3> <RNAcentral/bed> <genomeFASTA> <mature_miR_fasta> <related_sps_mature_fasta> <your_sps_precursor_fasta> <species 3 letter code>\n"
	exit 1;
elif [[ ${RNA_bed: -4} != ".bed" ]]; then
	echo -e "ERROR: RNACentral annotation is not a valid bed file format !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
else
	echo -e "#miRBase annotation:" "$miR_gff3"
	echo -e "#RNACentral annotation:" "$RNA_bed"
fi

echo -e "#Annotation file testing done, analysis proceeded\n"
# # echo -e	"#############################"
# # echo -e	"# # # ncRNA_FILTERATION # # #"
# # echo -e	"#############################"

if [[ -f $miR_anno ]] && [[ -f $RNAcentral ]]
then
	#echo -e "#ncRNA_FILTERATION started..\n"
	miR_dir=$(dirname "${miR_anno}") ## extract path only
	# # RNA_dir=$(dirname "${RNAcentral}") ## extract path only
	cat $miR_anno | grep -w 'miRNA' | grep 'chr*' | grep -v '#' | awk -v OFS='\t' '{ print $1,$4,$5,$3}' > $miR_dir/$sps.mir.gff3.bed
	cat $RNAcentral | grep 'chr*' > $miR_dir/$sps.rna.bed
	cat $miR_dir/$sps.rna.bed $miR_dir/$sps.mir.gff3.bed | awk -v OFS='\t' '{ print $1,$2,$2,$4}' | sort -k1,1 -k2,2n > $miR_dir/$sps.rna.mir.bed
	echo -e "#Filtering known annotations from bam.."
	for i in "${bam[@]}"; do
		# # filtration from miRBase genome co-ordinates
		intersectBed -abam ${i} -b $miR_dir/$sps.rna.mir.bed -wa -v > $input_map/mirdeep_input/${i##*/}.filtered.bam

		# # converting bam to sam for novel miRNA prediction
		##samtools view -h $input_map/mirdeep_input/${i##*/}.filtered.bam > $input_map/mirdeep_input/${i##*/}.sam
		# # discarding less than 18nt mapped reads and write only mapped reads
		samtools view -h -F 4 $input_map/mirdeep_input/${i##*/}.filtered.bam | awk '/^@/ || length($10) >= 18' > $input_map/mirdeep_input/${i##*/}.sam
		echo -e "Known RNA Transcript filtration has done for "${i##*/}"\n"
	done
else
	echo -e "Either miRBase or RNAcentral annotation files missing or Not correct file format,
please check input,analysis halted !!\n"
exit;
fi

# # declare an array variable
declare -a sam=($input_map/mirdeep_input/*.sam)

if [[ -f "./Scripts/bwa_sam_converter_edited.pl" ]]
then
	echo -e "#miRDeep2 inputs preparation started.."
	for j in "${sam[@]}"
	do
		filename="${j##*/}"
		bs=$(echo "$filename" | sed 's/\.sorted.bam.sam//')
		# # echo $bs
		# # convert sam into reads_collapsed and arf file format (inputs for miRDeep2)
		perl ./Scripts/bwa_sam_converter_edited.pl -i ${j} -c -o $input_map/mirdeep_input/${bs##*/}.reads_collapsed.fa -a $input_map/mirdeep_input/${bs##*/}.reads_collapsed_vs_genome.arf
		# # echo -e "${j##*/}"
		# # replacement by species
		sed -i "s/seq/$sps/g" $input_map/mirdeep_input/${bs##*/}.reads_collapsed.fa

		sed -i "s/seq/$sps/g" $input_map/mirdeep_input/${bs##*/}.reads_collapsed_vs_genome.arf

		# # discard reads less than 18 nt (if any)
		echo -e "#Discarding reads less than 18nt, miRDeep2 requirement.\n"
		perl ./Scripts/fastaparse.pl $input_map/mirdeep_input/${bs##*/}.reads_collapsed.fa -a 18 -b > $input_map/mirdeep_input/${bs##*/}.reads_collapsed_novel.fa 2> /dev/null
		echo -e "#miRDeep2 input has prepared for sample: "${bs##*/}.reads_collapsed.fa"\n"
	done
else
	echo -e "make sure you are in CBS-miRSeq directory,
miRDeep2 input preparation failed !!\n"
fi

echo -e "##################################"
echo -e "# # # Novel miRNA prediction # # #"
echo -e "##################################"

# declare an array variable
declare -a fa=($input_map/mirdeep_input/*novel.fa)
declare -a arf=($input_map/mirdeep_input/*.arf)

# length
readsLen=${#fa[@]}

for (( i=0; i<${readsLen}; i++ ));
do
	rep="${fa[$i]##*/}"
	repo=$(echo "$rep" | sed 's/\.reads_collapsed_novel.fa//')
	mkdir -p $output_dir/Novel_miR ## create dir name of each replicates (samples)
	mkdir -p $output_dir/Novel_miR/$repo
	cd $output_dir/Novel_miR/$repo ## must be go inside output dir dor mirdeep2 results
	echo -e "\n################ Begins:" `date` "############\n"
	echo -e "#Started miRDeep2 prediction for sample ::" "$repo"
	## -P = (5p and 3p) instead of previous ids
	## -c = disable randfold analysis
	## -g = -1 all precursors will be analyse
	## -u = hyper link of UCSC browser species that are supported
	miRDeep2.pl "${fa[$i]}" $input_map/mirdeep_input/genome_clean.fa "${arf[$i]}" $input_map/mirdeep_input/miR_clean.fa $input_map/mirdeep_input/related_miR_clean.fa $input_map/mirdeep_input/pre_mir_clean.fa -g -1 -t $species -u -c -v -P 2> $output_dir/Novel_miR/$repo/$repo.report.log
	echo -e "################ Finished:" `date` "##########""\n"
done


#################################################
## Remove all unused and temp dir after analysis
#################################################

# # # # # remove created annotation bed
rm -f $miR_dir/$sps.rna.bed
rm -f $miR_dir/$sps.mir.gff3.bed
rm -f $miR_dir/$sps.rna.mir.bed
# remove input dir
rm -rf $input_map/mirdeep_input
# remove temp dir of mirdeep2
results=$(find $output_dir/Novel_miR -name 'result*.csv'| wc -l)
if [[ $results -ne 0 ]]; then
	find $output_dir/Novel_miR -name 'dir_prepare_signature*' -type d -exec rm -r {} +
	find $output_dir/Novel_miR -name 'mirdeep_runs' -type d -exec rm -r {} +
	find $output_dir/Novel_miR -name 'expression_analyses' -type d -exec rm -r {} +
	find $output_dir/Novel_miR -name 'error_*' -exec rm -r {} +
else
	echo "Error, results not found: miRDeep2 has not predicted Novel miRs!
,Please check all inputs and run script again.\n"
fi

echo -e "\n_Novel miRNA Analysis finished": `date`

##================================= novel mirna parsing script wrapper ============##
echo -e "#=============Preparing a Novel matrix for further DE analysis===========#"

#################################### Analysis started #######################################

# # Step 1
# # header for novel matrix
# # find results csv from miRDeep2
declare -a  csv=($(find $output_dir/Novel_miR -name 'result_*.csv' | sort))

if [[ -z "${csv[0]}" ]]; then
    echo -e "csv does not exist\n"
	exit 1;
else
	for i in "${csv[@]}"; do
		awk '$13 !~ /-5p./{print}' $i > $i.tmp && mv $i.tmp $i ### remove undefined character created by miRDeep
		awk '$13 !~ /-3p./{print}' $i > $i.tmp && mv $i.tmp $i ### remove undefined character created by miRDeep
		# # header for novel matrix
		echo -n -e "${i##*/}\t" >> $output_dir/Novel_miR/novel.miRNA.counts.matrix.txt
		# # extract novel mirna, score and their read counts and then sort as bed
		##cat $i | sed -n '/provisional/,/^$/p' | awk '$2 >= 1' | sed '1d' | sed -e 's/:/\t/g;s/\../\t/g;s/ /\t/g' | awk '{OFS="\t"; print $19,$20,$21,$1,$2,$22,$7}' | sort -k1,1 -k2,2n > $output_dir/Novel_miR/${i##*/}.sorted.csv.bed
		cat $i | sed -n '/provisional/,/^$/p' | awk '$2 >= 1' | sed '1d' | sed -e 's/:/\t/g;s/\../\t/g;s/ /\t/g' | awk '{OFS="\t"; print $19,$20,$21,$1,$2,$22,$7}' > $output_dir/Novel_miR/${i##*/}.csv.bed
		bedtools sort -i $output_dir/Novel_miR/${i##*/}.csv.bed > $output_dir/Novel_miR/${i##*/}.sorted.csv.bed
		rm -rf $output_dir/Novel_miR/${i##*/}.csv.bed
		# # # # overlapped merge (output=chr start end counts)
		bedtools merge -i $output_dir/Novel_miR/${i##*/}.sorted.csv.bed -c 7 -o distinct > $output_dir/Novel_miR/${i##*/}.final.novel.bed.tmp
		# # # # give weight if merge reports more than one counts
		cat $output_dir/Novel_miR/${i##*/}.final.novel.bed.tmp | awk -v OFS="\t" '{n=split($NF, a, ","); for (i=1;i<=n;i++) s+=a[i]; $NF=s; s=0}1' > $output_dir/Novel_miR/${i##*/}.wt

	done
	# # concatenate into one file to create ref file
	##paste --delimiter=\\n --serial $output_dir/Novel_miR/*.bed | awk '{OFS="\t"; print $1,$2,$3}' | sort -k1,1 -k2,2n > $output_dir/Novel_miR/ref.file.sort
	paste --delimiter=\\n --serial $output_dir/Novel_miR/*.bed | awk '{OFS="\t"; print $1,$2,$3}' > $output_dir/Novel_miR/ref.file
	bedtools sort -i $output_dir/Novel_miR/ref.file > $output_dir/Novel_miR/ref.file.sort
	#rm -rf $output_dir/Novel_miR/ref.file
	# # merge to create ref file (create new ordinates because of repetition )
	bedtools merge -i $output_dir/Novel_miR/ref.file.sort > $output_dir/Novel_miR/ref.file.sort.refer
	echo "job 1 done"
fi

## To check if ref file has created
if [[ ! -e $output_dir/Novel_miR/ref.file.sort.refer ]]; then
	echo -e "*.ref does not exist; something was wrong in the processing of step one."
	echo -e "#Abort further processing !!\n"
	exit 1;
fi

# # # formatting header
sed 's/$/\n/g' $output_dir/Novel_miR/novel.miRNA.counts.matrix.txt | awk -v OFS="\t" '{print "#Chr","Start","End",$0}' > $output_dir/Novel_miR/tmp && mv $output_dir/Novel_miR/tmp $output_dir/Novel_miR/novel.miRNA.counts.matrix.txt 2> /dev/null

# Step 2
# match each sample with ref.file for counts

tmp=($output_dir/Novel_miR/*.wt)

if [[ -f "${tmp[0]}" ]]; then

for ref in $output_dir/Novel_miR/ref.file.sort.refer
do
   for file in $(ls -1 -t -v $output_dir/Novel_miR/*.wt)
   do
       if [[ $ref != $file ]]; then
	    intersectBed -a $ref -b $file -f 0.9 -r -wao | awk '$8>1{print $1,$2,$3,$7}' > $output_dir/Novel_miR/${ref##*/}.${file##*/}.common
       fi
   done
done
else
	echo -e "tmp files not founds..check errors !!!\n"
exit 1;
fi

## To check if *.common exist
cmm=($output_dir/Novel_miR/*.common)

if [[ -e "${cmm[0]}" ]]; then
##RENAMING
	for j in `ls $output_dir/Novel_miR/*.common`; do mv "$j" "${j/ref.file.sort.refer./}"; done
	echo -e "Step 2 done.\n"
else
	echo -e "common features files not found !!\n"
	exit 1;
fi

# #Step 3
# #make a expr matrix for diff expr analysis by perl script

# #Calling perl script

novel_parser="./Scripts/CBS-miRSeq.make_novel.miRNA.matrix_v1.0.pl"

if [[ -f $novel_parser ]]; then

	perl $novel_parser $output_dir/Novel_miR | sort -k1,1 -k2,2n >> $output_dir/Novel_miR/novel.miRNA.counts.matrix.txt
else
	echo -e "please make sure that your are in same directory of CBS-miRSeq directory.\n"
fi


# ##sleep 5
# delete tmp file
rm -f $output_dir/Novel_miR/*.tmp $output_dir/Novel_miR/*.common $output_dir/Novel_miR/*.wt $output_dir/Novel_miR/*.refer $output_dir/Novel_miR/*.ref

rm -f $output_dir/Novel_miR/*.sort $output_dir/Novel_miR/*.bed.tmp $output_dir/Novel_miR/ref.file


# # make header for bed files

declare -a bed=($output_dir/Novel_miR/*.bed)

if [[ -f "${bed[0]}" ]]; then

	for k in "${bed[@]}"; do

		sed -i -e "1i #Chr\tStart\tEnd\tProviID\tmiRDeep2Score\tStrand\tCounts" $k

	done
else
	echo -e "bed files are not found for writing header\n"
fi

echo -e "\n_Novel miRNA Parsing finished, created a Matrix": `date`
echo -e "Novel miRNA prediction has finished, go to your output.dir.path to check your results.\n"
clear
) 2>&1) | tee -a module3a.log

if [[ -e $logs/module3a.log ]]; then

	rm -f $logs/module3a.log

fi

mv -f module3a.log logs

echo -e "\n#log file also has created into current directory of CBS-miRSeq: logs/module3a.log\n"


exit $?

