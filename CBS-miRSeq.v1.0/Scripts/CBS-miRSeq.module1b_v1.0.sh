#!/bin/bash -xv

((
#=====================================CBS-miRSeq Module 1 analysis steps ================================================================================================================================
# Description: analysis pipeline for SOLiD reads (csfasta & qual) and Illumina reads (fastq)
# USAGE: bash ./CBS-miRSeq.module1b_v1.0.sh <input_reads_path> <output_path> <Index_dir> <miRBase annotation/sps.GFF3> <Ensembl annotation/sps.GTF> <species 3 letter code> <reads_color=yes or no> <genome_release(two digit integer)> <my_adaptor/base-space_seq or color-space>
# input_reads_dir_path = directory of csfasta and qual files
# output_dir_path = directory where results should be written
# bowtie_index_path = directory of genome index built by bowtie
# miRBase annotation/GFF3 = miRBase miRNA genome co-ordinates file in gff3
# Ensembl annotation/GTF = ensembl annotation in gff (genome co-ordinates)
# species 3 letter code = reference species 3 letter code (i.e. dre/mmu/hsa/rno etc)
# reads_color = "yes" if you have color-space reads from SOLiD (*.csfasta and qual) or otherwise "no" for illumina fastq reads
# Ens_release="84" ### ensembl release (two digit integer)
# my_adaptor="sequence of adaptor" ### a small RNA adaptor should be provide by the user
# Date : 15 May 2015
# version: 1.0
# Authors: Kesharwani RK email: bioinforupesh2009.au@gmail.com
#===========================================================================================================================================================================================================

#clear
echo -e '\0033\0143'
print_softInfo ()
{
echo -e "#^^^^^CBS-miRSeq Module 1 analysis steps of SOLiD reads (csfasta & qual) and Illumina (fastq)^^^^^#"
echo -e "#				~ ~ ~ CBS-miRSeq Module 1b v1.0 ~ ~ ~				  #"
echo -e "#		Analysis    : QC,Trimming,Mapping,Quantification				  #"
echo -e "#		Result: Summary statistics, mapping plots and Counts Matrix			  #"
echo -e "#		Requested Citation: Kesharwani RK et al.(2019)			  #"
echo -e "#		Author			  : bioinforupesh2009.au@gmail.com					  #"
echo -e "#		Copyright (c) 2019 Kesharwani RK				 #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
}

################################################################
# # exp: # # Note: use without quotes while using terminal
# # input_reads_dir="/data/..../reads" ## directory
# # output_dir="/data/...../results_dir" ## directory
# # bowtie_index="/data/.../Index" ## directory
# # miR_anno="/data/.../hsa.gff3" ## file
# # ensembl_anno="/data/.../Homo_sapiens.GRCh38.78.gtf" ## file
# # sps="hsa" ## string
# # reads_color="no" ## yes or no ## string
# # release="84 " ### ensembl release (two digit integer)
# # my_adaptor="sequence of adaptor" ### small RNA adaptor should be provide by user
###############################################################

print_softInfo

echo -e "\n"
echo -e "\t\t\t\t~ ~ ~ ~ ~ CBS-miRSeq.module1b_v1.0 ~ ~ ~ ~ ~"
echo -e "\t\t\t~ ~ ~ ~ ~ QC,Trimming,Mapping,Quantification & Summary ~ ~ ~ ~ ~"
echo -e "\t\t\t\tAnalysis date:" `date`
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

## Inputs: command line arguments
input_reads_dir="$1"
output_dir="$2"
bowtie_index="$3"
miR_anno="$4"
ensembl_anno="$5"
sps="$6"
reads_color="$7" ## yes or no
release="$8" ### ensembl release (two digit integer)
my_adaptor="$9" ### small RNA adaptor should be provide by user
NumberOfSamples="${10}"

## sps input checker
ref_cnt=$(echo -n "$sps" | wc -m)

print_USAGE()
{
	echo -e "USAGE: 
	bash ./CBS-miRSeq.module1b_v1.0.sh <input_reads_path> <output_path> <Index_dir> <miRBase annotation/sps.GFF3> <Ensembl annotation/sps.GTF> <species 3 letter code> <reads_color=yes or no> <genome_release(two digit integer)> <my_adaptor/base-space_seq or color-space> <NumberOfSamples>\n"
	echo -e "EXAMPLE: 
	bash ./CBS-miRSeq.module1b_v1.0.sh /.../reads /.../Results /.../Index /.../hsa.gff3 /.../Homo_sapiens.GRCh38.84.gtf hsa no 84 ATCTCGTATGCCGTCTTCTGCTTG 10"
	echo -e "\n"
}

##checking input arguments
if [[ $# -ne 10 ]]; then
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit;
elif [[ ! -d "$input_reads_dir" ]]; then
    echo -e "\n... Input directory "$input_reads_dir" is not a valid directory !!"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ $ref_cnt -ne 3 ]]; then
	my_color_text "#Please provide a 3 letter code of your reference species !!" red
	print_USAGE
	exit 1;
else
	echo -e "job running...\n"
	echo -e "File Testing....\n"
fi

#### checking bowtie index input
index=($bowtie_index/*.ebwt)
if [[ -e "${index[0]}" ]]; then
	echo -e "All Inputs are provided,analysis preceded.....\n"
	echo -e "#bowtie index exist"
else
	echo -e "*.ebwt(bowtie index) does not found!!
please provide a valid directory of input reads....
OR pre-built bowtie index of your species genome using CBS-miRSeq.module1a_v1.0.sh script,
before running this module !!\n"
exit 1;
fi

# # checking color space or base space files existence in "input_reads_dir"

csfasta=($input_reads_dir/*.csfasta)
qual=($input_reads_dir/*.qual)
fq=($input_reads_dir/*.fastq)

if [[ $reads_color == "yes" ]]
then
	if [[ -f "${csfasta[0]}" ]] && [[ -f "${qual[0]}" ]]; then
		echo -e "#color-space inputs are found."
	fi
elif [[ $reads_color == "no" ]]; then
	if [[ -f "${fq[0]}" ]]; then
		echo -e "#base-space inputs are found."
	fi
else
	echo -e "*.csfasta and *.qual or fastq files not found!!
please provide valid directory of input reads....\n"
exit 1;
fi

# # miRBase gff3 and ensembl annotations gtf file existence validation

miR_gff3=${miR_anno##*/}
RNA_gtf=${ensembl_anno##*/}
# miR_gff3 validation
if [[ ${miR_gff3: -5} != ".gff3" ]]; then
	echo -e "ERROR: miRBase annotation is not a valid gff3 file format !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ ${RNA_gtf: -4} != ".gtf" ]]; then
	echo -e "ERROR: Ensembl annotation is not a valid bed file format !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
else
	echo -e "#miR annotation:" "$miR_gff3"
	echo -e "#ensembl annotation:" "$RNA_gtf"
fi

# # ensembl release validation
rel=$(echo -n "$release" | wc -m) ## counts numeric character which should be 2
if [[ "$rel" -ne 2 ]]; then
	echo -e "\n#ERROR: provided "$release" is not a correct integers,
please provide 2 digit integers only, make sure ensembl has this "$release""
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
exit 1;
else
	echo -e "#Species:" $sps
	echo -e "#Ensembl release:" $release"\n"
fi

echo -e "File testing done, analysis being analysed......\n"

# # # # # Make a directory for converted reads and their QC
mkdir -p $input_reads_dir/reads
## analysis directory
mkdir -p $output_dir/{before_qc,after_qc,clean_reads,discarded_reads,trim_reads,log_files,reads_mapped,visualization_summary,quantification,DE}

# # ## make user's reads back up and used this for analysis
### checking color-space and base-space files existence in "reads" directory
csfasta=($input_reads_dir/reads/*.csfasta)
qual=($input_reads_dir/reads/*.qual)
fq=($input_reads_dir/reads/*.fastq)

if [[ $reads_color == "yes" ]]
then
	if [[ ! -f "${csfasta[0]}" ]] && [[ ! -f "${qual[0]}" ]]
	then
		echo -e "cleaning invisible carriage and making a backup of color-space reads are in progress, please wait..."
		## remove blank and comments lines if there is:
		# # ## declare an array variable for input
		declare -a  in_reads=($input_reads_dir/*.csfasta)
		for i in "${in_reads[@]}"; do sed -i '/^$/d' $i | grep -v -e'^#'; done
		declare -a  inq_reads=($input_reads_dir/*.qual)
		for j in "${inq_reads[@]}"; do sed -i '/^$/d' $j | grep -v -e'^#'; done
		cp $input_reads_dir/*.{csfasta,qual} $input_reads_dir/reads
		echo -e "back up done..\n"
	fi
elif [[ $reads_color == "no" ]]
then
	if [[ ! -f "${fq[0]}" ]]
	then
		echo -e "cleaning invisible carriage and making a backup of base-space reads are in progress, please wait..."
		##remove blank and comments lines if there is:
		# # ## declare an array variable for input
		declare -a  in_reads=($input_reads_dir/*.fastq)
		for i in "${in_reads[@]}"; do sed -i '/^$/d' $i; done
		cp $input_reads_dir/*.fastq $input_reads_dir/reads
		echo -e "back up done..\n"
	fi
else
	echo "can not make back-up, analysis halted !! \n"
exit 1;
fi

echo -e "\t\t\t#^^^^^CBS-miRSeq Module 1 analysis started^^^^^#\n"


echo -e "####################################"
echo -e "## QC analysis before filtering   ##"
echo -e "####################################"

echo -e "\nStarted...\n"

if [[ $reads_color == "yes" ]]
then
	# #read all file name into an array
	cs=($(find $input_reads_dir/reads/ -name *.csfasta | sort))
	qv=($(find $input_reads_dir/reads/ -name *.qual | sort))
	# #get length of an array
	cLen=${#cs[@]}
	# #use for loop read all file names
	for (( i=0; i<${cLen}; i++ ));
	do
		# # basename
		base="${cs[$i]##*/}"
		##echo "$base"
		##echo "${cs[$i]}"
		echo -e "Color-space Conversion for before QC...\n"
		# # echo -e "${cs[$i]}" "${qv[$i]}"
		solid2fastq -o > /dev/null 2>&1 $output_dir/before_qc/${cs[$i]##*/} "${cs[$i]}" "${qv[$i]}" > /dev/null
		echo -e "Before QC for ${cs[$i]##*/}\n"
		fastqc -o $output_dir/before_qc/ --noextract -q $output_dir/before_qc/${cs[$i]##*/}.fastq
		echo -e "done both\n."
	done
elif [[ $reads_color == "no" ]]
then
	fqr=($input_reads_dir/reads/*.fastq)
	for j in "${fqr[@]}"; do
	echo -e "Before QC for base space started for ${j##*/}...\n"
	fastqc -o $output_dir/before_qc/ --noextract -q $j
	echo -e "base space QC done for ${j##*/}\n"
	done
else
	echo "can not perform filtration and before QC, analysis halted !! \n"
exit 1;
fi
echo -e "\nDone.\n"

echo -e "#################################"
echo -e "## Read quality filtering (QC) ##"
echo -e "#################################"
echo -e "\nStarted...\n"

filt_script="SOLiD_preprocess_filter_v2.pl"
if [[ $reads_color == "yes" ]]
then
	if [[ -f "$filt_script" ]]; then
		# #read all file name into an array
		cs=($(find $input_reads_dir/reads/ -name *.csfasta | sort))
		qv=($(find $input_reads_dir/reads/ -name *.qual | sort))
		# #get length of an array
		cLen=${#cs[@]}
		# #use for loop read all file names
		for (( i=0; i<${cLen}; i++ ));
		do
			# # basename
			base="${cs[$i]##*/}"
			# #echo "$base"
			# #echo "${base%%.*}"
			# #echo "${cs[$i]}"
			echo -e "color reads filtration begin for ${cs[$i]##*/} ${qv[$i]##*/}...\n"
			perl $filt_script -f "${cs[$i]}" -g "${qv[$i]}" -x off -y on -e 3 -d 9 -n on -v on -o $output_dir/clean_reads/"${cs[$i]##*/}" > /dev/null
			echo -e "color reads filtration done.\n"
		done
	else
		echo "you are not in CBS-miRSeq scripts directory, analysis halted !!\n"
	exit 1;
	fi
elif [[ $reads_color == "no" ]]; then
	fqr=($input_reads_dir/reads/*.fastq)
	for j in "${fqr[@]}"; do
		name="${j##*/}"
		##echo "${name%.*}"
		echo -e "base-space filtration begins: For ${j##*/}...\n"
		if grep -q '^Z' "$j"; then
			echo -e "it is a phred64 quality reads !!"
			fastq_quality_filter -q 20 -p 50 -Q 64 -v -i $j -o $output_dir/clean_reads/"${name%.*}".cln.fastq > $output_dir/log_files/"${name%.*}".cln.log
		else
			echo -e "it is a phred33 quality reads !!"
			fastq_quality_filter -q 20 -p 50 -Q 33 -v -i $j -o $output_dir/clean_reads/"${name%.*}".cln.fastq > $output_dir/log_files/"${name%.*}".cln.log
			echo -e "QC done: For ${j##*/}\n"
		fi
	done
	# # # filtered log join together
	cat $output_dir/log_files/*.cln.log > $output_dir/log_files/fastq_reads.log
else
	echo -e "reads not founds to be filter..!! \n"
exit 1;
fi
echo -e "\nDone.\n"

# # # moving discarded reads
mv $output_dir/clean_reads/*_U_F3* $output_dir/discarded_reads/ 2> /dev/null

echo -e "###################################"
echo -e "## Trimming (3' adaptor removal) ##"
echo -e "###################################"
echo -e "\nStarted...\n"

csfasta=($output_dir/clean_reads/*.csfasta)
qual=($output_dir/clean_reads/*.qual)

if [[ $reads_color == "yes" ]]
then
	if [[ -f "${csfasta[0]}" ]] && [[ -f "${qual[0]}" ]]
	then
		echo -e "QC reads found...Trimming proceeded.. \n"
		# #read all file name into an array
		cs_cln=($(find $output_dir/clean_reads -name *.csfasta | sort))
		qv_cln=($(find $output_dir/clean_reads -name *.qual | sort))

		# #get length of an array
		cLen=${#cs_cln[@]}

		# #use for loop read all file names
		for (( i=0; i<${cLen}; i++ ));
		do
			# # basename
			base="${cs_cln[$i]##*/}"
			##echo "${base%%.*}"
			# # print file
			#echo "${cs_cln[$i]}" "${qv_cln[$i]}"
			cs_adaptor="330201030313112312" #default adaptor
			if [[ "$cs_adaptor" != "$my_adaptor" ]]; then
				echo -e "Trimming begin for color reads "${cs_cln[$i]##*/}" "${qv_cln[$i]##*/}"...\n"
				echo -e "#user's adaptor is using to trim the reads."
				cutadapt -z -c -e 0.2 -m 15 -a "$my_adaptor" "${cs_cln[$i]}" "${qv_cln[$i]}" -o $output_dir/trim_reads/${cs_cln[$i]##*/}.fastq > $output_dir/log_files/${base%%.*}.trim.tmp.log
			else
				echo -e "#default and provided adaptor is the same; Trimming proceeded.."
				echo -e "Trimming begin for color reads "${cs_cln[$i]##*/}" "${qv_cln[$i]##*/}"...\n"
				cutadapt -z -c -e 0.2 -m 15 -a "$cs_adaptor" "${cs_cln[$i]}" "${qv_cln[$i]}" -o $output_dir/trim_reads/${cs_cln[$i]##*/}.fastq > $output_dir/log_files/${base%%.*}.trim.tmp.log
			fi
			echo -e "FastQC after filter started...\n"
			# #QC analysis after trimmed reads
			fastqc -o $output_dir/after_qc --noextract -q $output_dir/trim_reads/${cs_cln[$i]##*/}.fastq
			echo -e "Trimming and QC done successfully.. for "${cs_cln[$i]##*/}"\n"
		done
	else
		echo -e "color clean reads not found, Trimming and QC halted !! \n"
	fi
elif [[ $reads_color == "no" ]]; then
	fqcln=($output_dir/clean_reads/*.fastq)
	for cln in "${fqcln[@]}"; do
		name="${cln##*/}"
		# #echo "${name%.*}"
		echo -e "base-space trimming begins for ${cln##*/}...\n"
		# #echo "${name%.*}"
		## v1.5 Small RNA 3' adaptor (ATCTCGTATGCCGTCTTCTGCTTG)/Illumina Single End adaptor 1
		## http://supportres.illumina.com/documents/documentation/chemistry_documentation/experiment-design/illumina-customer-sequence-letter.pdf
		## hsa miRBase 21: min length of mature miR = 16 nt, max = 28 nt
		## Hence, by default min length of reads are = 15nt and max = no limit
		bs_adaptor="ATCTCGTATGCCGTCTTCTGCTTG" #default adaptor
		if [[ "$bs_adaptor" != "$my_adaptor" ]]; then
			echo -e "#user's adaptor is using to trim the reads."
			cutadapt -e 0.2 -m 15 -a "$my_adaptor" "$cln" -o $output_dir/trim_reads/"${name%.*}".fastq > $output_dir/log_files/"${name%.*}".trim.tmp.log
		else
			echo -e "#default and provided adaptor is the same; Trimming proceeded.."
			cutadapt -e 0.2 -m 15 -a "$bs_adaptor" "$cln" -o $output_dir/trim_reads/"${name%.*}".fastq > $output_dir/log_files/"${name%.*}".trim.tmp.log
		fi
		# #QC analysis after trimmed reads
		fastqc -o $output_dir/after_qc --noextract -q $output_dir/trim_reads/"${name%.*}".fastq
		echo -e "base space trimming and after QC done for ${cln##*/}\n"
	done
else
	echo -e "reads not founds,Trimming and 2nd QC halted..!! \n"
exit 1;
fi
echo -e "\nDone.\n"

# # ## join all adaptor logs into one file and move it for summary processing
cat $output_dir/log_files/*.trim.tmp.log > $output_dir/log_files/all.trim.log 2> /dev/null

echo -e "#############"
echo -e "## Mapping ##"
echo -e "#############"
echo -e "\nStarted...\n"
# # trimmed reads check
t_reads=($output_dir/trim_reads/*.fastq)
m_reads=($output_dir/reads_mapped/*.sam)

if [[ $reads_color == "yes" ]]
then
	if [[ -f "${t_reads[0]}" ]]
	then
		echo -e "Trimmed reads found...mapping proceeded.. \n"
		if [[ ! -f "${m_reads[0]}" ]]
		then
			fqtrm=($output_dir/trim_reads/*.fastq)
			for trim_reads in "${fqtrm[@]}"; do
				### basename
				rbname=$(basename "$trim_reads" .csfasta.fastq)
				echo -e "color-space Mapping begins..., for ${trim_reads##*/}\n"
				bowtie --sam --best --col-keepends --chunkmbs 256 -C -q -p 4 -n 0 -l 15 $bowtie_index/$sps.genome_v"$release"C ${trim_reads} $output_dir/reads_mapped/${rbname}.sam 2> $output_dir/log_files/${rbname}.map.tmp.log
				### converting sam directly to a sorted bam
				samtools view -bS $output_dir/reads_mapped/${rbname}.sam | samtools sort - > /dev/null $output_dir/reads_mapped/${rbname}.sorted
				## creating bam index of aligned reads (*.bai)
				samtools index $output_dir/reads_mapped/${rbname}.sorted.bam > /dev/null
				### Mapping statistics and plots as html output
				samstat $output_dir/reads_mapped/${rbname}.sam 2> /dev/null
				echo -e "color space Mapping done for ${trim_reads##*/}\n"
			done
		fi
	else
		echo -e "color space trimmed reads not found..\n"
	fi
elif [[ $reads_color == "no" ]]; then
	fqtrm=($output_dir/trim_reads/*.cln.fastq)
	for trim_reads in "${fqtrm[@]}"; do
		### basename
		rbname=$(basename "$trim_reads" .cln.fastq)
		echo -e "base-space mapping begins: For ${trim_reads##*/}\n"
		bowtie --sam --best --chunkmbs 256 -q -p 4 -n 0 -l 15 $bowtie_index/$sps.genome_v"$release"B ${trim_reads} $output_dir/reads_mapped/${rbname}.sam 2> $output_dir/log_files/${rbname}.map.tmp.log
		### converting sam directly to a sorted bam
		samtools view -bS $output_dir/reads_mapped/${rbname}.sam | samtools sort - > /dev/null $output_dir/reads_mapped/${rbname}.sorted
		## creating bam index of aligned reads (*.bai)
		samtools index $output_dir/reads_mapped/${rbname}.sorted.bam > /dev/null
		### Mapping statistics and plots as html output
		samstat $output_dir/reads_mapped/${rbname}.sam 2> /dev/null
		echo -e "mapping done: For ${trim_reads##*/}\n"
	done
else
	echo -e "Trim fastq reads not found, alignment and other analysis halted !!\n"
exit 1;
fi
echo -e "\nDone.\n"
### move samstats files into visualization directory
mv -f $output_dir/reads_mapped/*.html $output_dir/visualization_summary 2> /dev/null

##join all mapping logs into one file and move it for summary processing
cat $output_dir/log_files/*.map.tmp.log > $output_dir/log_files/all.map.log 2> /dev/null

echo -e "#####################"
echo -e "## Quantification ###"
echo -e "#####################"
echo -e "\nStarted...\n"

touch "$output_dir/quantification/$sps.ensembl.report.txt" 2> /dev/null
touch "$output_dir/quantification/$sps.miRBase.report.txt" 2> /dev/null

## Run featuresCounts

map=`find $output_dir/reads_mapped/ -name *.sam | sort`

declare -a map=($output_dir/reads_mapped/*.sam)

if [[ -f "$miR_anno" ]] && [[ -f "$ensembl_anno" ]];
then
	echo -e "Quantification started, please wait...\n"
	# # echo -e "$map" | xargs -I % basename "%" ".sam"
	# # Counts gene features from ensembl for classification
	featureCounts -T 2 -t gene -g gene_biotype -a $ensembl_anno "${map[@]}" -o $output_dir/quantification/$sps.ensembl.gene.class.tmp 2> /dev/null
	# ## Counts gene features from miRBase for Diff Expression analysis
	featureCounts -T 2 -t miRNA -g Name -a $miR_anno "${map[@]}" -o $output_dir/quantification/$sps.miR.counts.tmp 2> /dev/null
	echo -e "quantification done !!\n"
else
	echo -e "Reference annotation files are not found, quantification skipped..\n"
exit 1;
fi

##extract the filenames in order as featurescounts done
if [[ -f $output_dir/quantification/$sps.miR.counts.tmp.summary ]]; then
	plusoneSample=$(($NumberOfSamples+1))
	for i in $(seq 2 $plusoneSample); do 
		pat=($(cat $output_dir/quantification/$sps.miR.counts.tmp.summary | head -1 | cut -f$i))
		basename $pat .sam | tr "\n" "\t" >> $output_dir/quantification/$sps.miRBase.report.txt
		basename $pat .sam | tr "\n" "\t" >> $output_dir/quantification/$sps.ensembl.report.txt
	done
else
	echo -e "Output of features summary is empty.!!\n"
	exit 1;
fi

## creating miR counts report
# # formatting header
sed 's/$/\n/g' $output_dir/quantification/$sps.miRBase.report.txt | awk -v OFS="\t" '{print "miRNAs",$0}' > $output_dir/quantification/tmpM && mv $output_dir/quantification/tmpM $output_dir/quantification/$sps.miRBase.report.txt 2> /dev/null
## extract counts from count table generated by featuresCounts
if [[ -f $output_dir/quantification/$sps.miR.counts.tmp ]]; then
	# # Formatting
	sed '1,2d' $output_dir/quantification/$sps.miR.counts.tmp | awk '{OFS="\t";$2=$3=$4=$5=$6="";$0=$0;$1=$1}1' > $output_dir/quantification/$miR.counts.txt
	# # header
	cat $output_dir/quantification/$miR.counts.txt >>  $output_dir/quantification/$sps.miRBase.report.txt
	## clear temp file
	rm $output_dir/quantification/$miR.counts.txt
else
	echo -e "Output of miR Counts are empty.!!\n"
	exit 1;
fi

## creating ensembl biotypes report
## Formatting
sed 's/$/\n/g' $output_dir/quantification/$sps.ensembl.report.txt | awk -v OFS="\t" '{print "gene",$0}' > $output_dir/quantification/tmpE && mv $output_dir/quantification/tmpE $output_dir/quantification/$sps.ensembl.report.txt 2> /dev/null
## extract counts from count table generated by featuresCounts
if [[ -f $output_dir/quantification/$sps.ensembl.gene.class.tmp ]]; then
	# # Formatting
	cat $output_dir/quantification/$sps.ensembl.gene.class.tmp | sed '1,2d' | awk '{OFS="\t";$2=$3=$4=$5=$6="";$0=$0;$1=$1}1' | grep -w 'lincRNA\|miRNA\|misc_RNA\|Mt_rRNA\|Mt_tRNA\|snRNA\|snoRNA\|protein_coding\|rRNA\|processed_transcript\|pseudogene\|sense_overlapping' | sort -k 1 > $output_dir/quantification/$sps.ensembl.class.txt
	# # header
	cat $output_dir/quantification/$sps.ensembl.class.txt >> $output_dir/quantification/$sps.ensembl.report.txt
	## clear temp file
	rm $output_dir/quantification/$sps.ensembl.class.txt
else
	echo -e "Output of Ensembl Biotypes Counts are empty.!!\n"
	exit 1;
fi

echo -e "\nDone...\n"

echo -e "#####################################"
echo -e "## summary statistic of each step ###"
echo -e "#####################################"
echo -e "\nStarted...\n"

if [[ $reads_color == "yes" ]]
then
	# # Raw reads QC stats
	touch $output_dir/log_files/raw.cnts.txt
	rw_read=($output_dir/before_qc/*.fastq)
	# # raw reads counts
	if [[ -s "${rw_read[0]}" ]]
	then
		echo -e "raw reads counting begin..\n"
		grep -c -H '^T' $output_dir/before_qc/*.fastq | xargs -L 1 basename | sed 's/\:/\t/g;s/.fastq//g' >> $output_dir/log_files/raw.cnts.txt
		echo -e "raw reads counts done..\n"
	else
		echo -e "Raw fastq converted reads not founds !! raw reads counting skipped.\n"
	exit 1;
	fi

	touch $output_dir/log_files/dis.cnts.txt
	ds_read=($output_dir/discarded_reads/*.csfasta)

	# # discarded reads counts
	if [ -s "${ds_read[0]}" ]
	then
		echo -e "discarded reads counting begin..\n"
		grep -c -H '>' $output_dir/discarded_reads/*.csfasta | xargs -L 1 basename | sed 's/\:/\t/g' | cut -f2 >> $output_dir/log_files/dis.cnts.txt
		echo -e "discarded reads counts done..\n"
	else
		echo -e "filtered reads reads not founds !! counting skipped.\n"
		exit 1;
	fi
elif [[ $reads_color == "no" ]]
then
	# # Raw reads QC stats
	touch $output_dir/log_files/raw.cnts.txt
	echo -e "raw reads counting begin..\n"
	for j in $input_reads_dir/reads/*.fastq; do echo -en "${j##*/}\t"; expr $(cat $j | wc -l) / 4; done >> $output_dir/log_files/raw.cnts.txt
	echo -e "raw reads counts done..\n"

	# # discarded reads counts
	touch $output_dir/log_files/dis.cnts.txt
	grep -w "discarded" $output_dir/log_files/fastq_reads.log | awk '{print $2}' >> $output_dir/log_files/dis.cnts.txt
	echo -e "discarded reads counts done..\n"
fi

# # Trimmed statistics

if [[ -s "$output_dir/log_files/all.trim.log" ]]
then
	echo -e "Trim counting begin..\n"
	cat $output_dir/log_files/all.trim.log | grep 'Processed reads' | awk '{print $3}' > $output_dir/log_files/processed.txt ## afterQC clean reads
	cat $output_dir/log_files/all.trim.log | grep 'Too long reads' | awk 'OFS="";{print $4,$5}' | sed 's/%/%)/g' > $output_dir/log_files/too_long.txt
	cat $output_dir/log_files/all.trim.log | grep 'Too short reads' | awk 'OFS="";{print $4,$5}' | sed 's/%/%)/g' > $output_dir/log_files/too_short.txt
	cat $output_dir/log_files/all.trim.log | grep 'Trimmed read' | awk 'OFS="";{print $3,$4}' > $output_dir/log_files/trimmed.txt
	echo -e "Trim counts done..\n"
else
	echo -e "trimmed "all.trim.log" reads not founds !! counting skipped.\n"
exit 1;
fi

# # Mapping statistics
if [[ -s "$output_dir/log_files/all.map.log" ]]
then
	echo -e "Map counting begin..\n"
	cat $output_dir/log_files/all.map.log | grep 'alignment:' | awk 'OFS="";{print $9,$10}' > $output_dir/log_files/aligned.reads.txt
	cat $output_dir/log_files/all.map.log | grep 'reads processed:' | awk 'OFS="";{print $4}' > $output_dir/log_files/sent_to_map.txt
	echo -e "Map counts done..\n"
else
	echo -e "Map logs reads not founds !! counting skipped.\n"
exit 1;
fi

#  calculate percentage of filtered reads and paste altogether
paste $output_dir/log_files/{raw.cnts.txt,dis.cnts.txt,processed.txt,trimmed.txt,too_short.txt,too_long.txt,sent_to_map.txt,aligned.reads.txt} \
	| awk -v OFS="\t" '{for(i=3;i<5;i++)$i=sprintf("%d(%.2f%)", $i, ($i/$2)*100)}1' > $output_dir/log_files/summary.stats.xls

# header
echo -e "Sample\tTotal counts\tFiltered reads\tReads after QC\tTrimmed reads\tToo short reads(<15bp)\tToo long reads (>35bp)\tReads sent to mapping\tMapped reads" \
 | cat - $output_dir/log_files/summary.stats.xls > $output_dir/log_files/tmp && mv $output_dir/log_files/tmp $output_dir/log_files/summary.stats.xls

echo -e "Summary excel file created..\n"

# ## moving summary xls
mv $output_dir/log_files/summary.stats.xls $output_dir/visualization_summary 2> /dev/null
echo -e "\nDone.\n"


echo -e "#####################################################################"
echo -e "## creation of HTML page for expressed miRNAs with link to miRBase	##"
echo -e "#####################################################################"
echo -e "\nStarted...\n"

## To recognize file size,whether file contains miRNA counts

file="$output_dir/quantification/$sps.miRBase.report.txt"
minimumsize=50 ## empty file size (with header)
actualsize=$(wc -c <"$file")
if [[ $actualsize -ge $minimumsize ]]; then
		html_script="./Scripts/CBS-miRSeq.HTML_v1.R"
		if [[ -f "$html_script" ]]; then

			Rscript "$html_script" $output_dir/quantification/$sps.miRBase.report.txt $output_dir/quantification "$sps"

		else
			echo -e "CBS-miRSeq.HTML_v1.R not found !!\n probably you are not in CBS-miRSeq script directory"
		fi
	else
	echo -e "$sps.miRBase.report.txt" or counts "are empty, please check back to all input and relaunch this module 1b."
	exit 1;
fi



# # cleaning directories..
# # remove converted raw reads, it was also used for raw counting
rm -f $output_dir/before_qc/*.fastq 2> /dev/null
# # removing temp files
rm -rf $output_dir/log_files/*.{txt,trim.tmp.log,map.tmp.log,cln.log} 2> /dev/null
rm -rf $output_dir/quantification/*.{ensembl.class.txt,miR.counts.txt,tmp,summary} 2> /dev/null


echo -e "CBS-miRSeq Module 1 has completed,
go through your output directory and check results,
now you may proceed for CBS-miRSeq.module2a,2b\n"

echo -e "\nAnalysis finished": `date`


) 2>&1) | tee -a module1b.log

if [[ -s $logs/module1b.log ]]; then

	rm -f $logs/module1b.log

fi

mv -f module1b.log logs

echo -e "\n#log file also has created in the directory of CBS-miRSeq: logs/module1b.log\n"

exit $?


