#!/bin/bash

((
#========================Prerequisite step to run main CBS-miRSeq Modules=================================
# Description:	    building.index of species genome
# USAGE: bash ./CBS-miRSeq.module1a_v1.0.sh <genome_dir_path> <reference_species_3_letter_code> <index_color=yes or no> <Ensembl_release(two digit integer)> <GRCh(Ensembl_chromosomal_assembly)>
# genome_dir="/.../dir" #directory where genome should be download or you have it and further will use to build index
# sps="hsa" # reference species 3 letter code
# index_color="yes" #(for color space reads from SOLiD) or no (for base space reads from illumina)
# Ens_release="84" ### ensembl release (two digit integer)
# GRCh="38" ## ensembl chromosomal assembly as numeric value
# This script will download your species genome if you don not have it as well as will build bowtie index
# Please make sure you have enough space (>= 10 GB)
# Date : 15 May 2015
# version = v1.0
# Authors: Kesharwani RK email: bioinforupesh2009.au@gmail.com
#=========================================================================================================

#clear
echo -e '\0033\0143'
print_softInfo () {
echo ''
echo -e "#^^^^^^^^^^^^^^^Prerequisite step to run CBS-miRSeq Modules 1a^^^^^^^^^^^^^^^#"
echo -e "#		   Analysis:	  building.index of species genome	     #"
echo -e "#		    Results: $get Genome FASTA and build bowtie index	     #"
echo -e "#	 Requested Citation: Kesharwani RK et al.(2019)	     #"
echo -e "#		     Author:  Kesharwani RK; bioinforupesh2009.au@gmail.com	     #"
echo -e "#		Copyright (c) 2019 Kesharwani RK				 #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
echo ''
}
print_softInfo

echo -e "\n"
echo -e "\t\t\t	   ~ ~ ~ ~ ~ CBS-miRSeq.module1a_v1.0~ ~ ~ ~ ~"
echo -e "\t\t\t~ ~ ~ ~ ~building.index of species genome~ ~ ~ ~ ~"
echo -e "\t\t\t	  Analysis date:" `date`

#################################################################################
# # Example inputs: # # Note: use without quotes while using terminal
# # genome_dir="/....../genome"
# # sps="hsa" ## your analysis species
# # index_color="no" ## yes or no
# # Ens_release="84" ### ensembl release (two digit integer)
# # GRCh="38" ##ensembl chromosomal assembly as numeric value
# # example command for hsa (Human)
# # bash ./CBS-miRSeq.module1a_v1.0.sh /data/..../directory/ hsa no 84 38
#################################################################################

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


mkdir -p logs
#Inputs: command line arguments
genome_dir="$1" ## directory
sps="$2" # three letter code of reference species
index_color="$3" ## yes or no
release="$4" ### ensembl release (two digit integer)
GRCh="$5" ## chromosomal assembly as numeric value (such as 38 or 36.67)

## for checking length of sps
ref_cnt=$(echo -n "$sps" | wc -m)

print_USAGE()
{
echo -e "USAGE: bash ./CBS-miRSeq.module1a_v1.0.sh <genome_dir_path> <reference_species_3_letter_code> <index_color=yes or no> <Ensembl_release(two digit integer)> <GRCh(Ensembl_chromosomal_assembly)>\n"
my_color_text "EXAMPLE:" blue
my_color_text "bash ./CBS-miRSeq.module1a_v1.0.sh /....../genome/ hsa no 84 38" red
echo -e "\n"
}

# checking input arguments
if [[ $# -ne 5 ]]; then
	my_color_text "#please supply all inputs.." cyan
	echo -e "\n"
	print_USAGE
	echo -e "#==============================================================================================================#"
	my_color_text "NOTE: Be careful with the input of Ensembl release and chromosomal assembly" cyan
	my_color_text "visit for info: http://www.ensembl.org/info/data/ftp/index.html" cyan
	my_color_text "NOTE: genome_dir_path is a directory where genome must download (If does not exist) and build index for the alignment" blue
	echo -e "#==============================================================================================================#"
	exit;
elif [[ ! -d "$genome_dir" ]]; then
    echo -e "\n#Input directory "$genome_dir" \ndoes not exist !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ $ref_cnt -ne 3 ]]; then
	echo -e "\n#Please provide 3 letter code of your reference species !!\n"
	echo -e "\nError(2)... \nplease supply all inputs.. \n"
	print_USAGE
	exit 1;
elif [[ $index_color != "yes" ]] && [[ $index_color != "no" ]]; then
	echo -e "\n#Index flag not found !! should be only: yes or no, analysis halted.. !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
else
	echo -e "Inputs are fine !! testing proceeded..\n"
fi

#check OS (Unix/Linux or Mac)
os=`uname`;
if [[ "$os" = "Darwin" ]]; then
	# use curl as the download program
	get="curl -L -o"
else
	# use $get as the download program
	get="wget --no-check-certificate -O"
fi


# # ensembl release validation
rel=$(echo -n "$release" | wc -c) ## counts numeric character which should be 2
if [[ "$rel" -ne 2 ]]; then
	echo -e "\n#ERROR: provided "$release" is not a correct integers,
please provide 2 digit integers only, make sure ensembl has this "$release""
	echo -e "\nError(4)... \nplease supply all inputs.. \n"
	print_USAGE
	exit 1;
else
	echo -e "\n"
	echo -e "#Species: "$sps""
	echo -e "#Ensembl release: "$release""
	echo -e "#chromosomal assembly: "$GRCh""
	echo -e "#Building color index: "$index_color""
	echo -e "#Output directory: "$genome_dir""
fi

## checking Genome FASTA (*.fa) input if user already have it.
# Already present genome fasta should be from Ensembl ; bcz other source can be a problem due to header info
# http://www.ensembl.org/info/data/ftp/index.html
gfa=($genome_dir/*.fa)
gfasta=($genome_dir/*.fasta)
if [[ -f "$gfa" ]] || [[ -f "$gfasta" ]]; then
	echo -e "#Reference genome fasta is exist, skipping downloading."
	## renaming if genome already exist in directory
	mv $genome_dir/*.fa $genome_dir/$sps.genome.tmp.fa 2> /dev/null
	mv $genome_dir/*.fasta $genome_dir/$sps.genome.tmp.fa 2> /dev/null
	if grep -q '^>chr' "$genome_dir/$sps.genome.tmp.fa"
	then
		echo -e "#chr header is fine, formatting is not required\n"
		mv $genome_dir/$sps.genome.tmp.fa $genome_dir/$sps.genome_v$release.fa
	else
		echo -e "#chr header NOT fine, formatting header\n"
	## formatting header
		awk '{print $1}' $genome_dir/$sps.genome.tmp.fa > $genome_dir/$sps.genome_v$release.fa
		sed -i -e 's/^>/&chr/g' $genome_dir/$sps.genome_v$release.fa
	fi
	sleep 2
	rm -f $genome_dir/$sps.genome.tmp.fa 2> /dev/null
	if [[ $index_color == "yes" ]]; then
		echo -e "Attention:\tbuilding COLOR-SPACE INDEX for your csfasta reads..please wait !!\n"
		echo -e "This may take a while, sorry but we cant speedup the process for bowtie-build.\n"
		bowtie-build -C $genome_dir/$sps.genome_v$release.fa $genome_dir/$sps.genome_v"$release"C > /dev/null
		echo -e "color-space index built.\n"
		echo -e "\nAnalysis finished:" `date`
		echo -e "Index successfully built, now run CBS-miRseq module 1b for the further analysis..\n"
		mkdir -p $genome_dir/Index
		mv $genome_dir/*.{ebwt,fa,fasta} $genome_dir/Index 2>/dev/null
		cp $genome_dir/Index/*.{fa,fasta} $genome_dir 2>/dev/null
		exit 1;
	elif [[ $index_color == "no" ]]; then
		echo -e "Attention:\tbuilding BASE-SPACE INDEX for your fastq reads..please wait !!\n"
		echo -e "This may take a while, sorry but we cant speedup the process for bowtie-build.\n"
		bowtie-build $genome_dir/$sps.genome_v$release.fa $genome_dir/$sps.genome_v"$release"B > /dev/null
		echo -e "base-space index built.\n"
		echo -e "\nAnalysis finished:" `date`
		echo -e "Index successfully built, now run CBS-miRseq module 1b for the further analysis..\n"
		mkdir -p $genome_dir/Index
		mv $genome_dir/*.{ebwt,fa,fasta} $genome_dir/Index 2>/dev/null
		cp $genome_dir/Index/*.{fa,fasta} $genome_dir 2>/dev/null
		exit 1;
	else
		echo -e " index not build, please check inputs !!\n"
		exit 1;
	fi
else
	echo -e "#Reference genome FASTA(*.fa) file does not found, downloading proceeded.."
	echo -e "This may take a while, Depend of your system speed.\n"
fi

# checking supported sps for downloading corresponding genome from ensembl genome browser
if [[ $sps == "hsa" ]] || [[ $sps == "dre" ]] || [[ $sps == "mmu" ]] || [[ $sps == rno ]]; then
	echo -e "\n"
	echo -e "job running...\n"
else
	echo -e "Currently supported reference Species to fetch..\n"
	echo -e "hsa	Human
dre	Zebrafish
mmu	Mouse
rno	Rat\n"
    echo -e "Given <"$sps"> species does not supported to fetch,
please download it manually and RUN this script again to build bowtie index.!!
Source to download genome from Ensembl(only Chr and MT; mandatory and preferred by the CBS-miRSeq pipeline):
#Ensembl Release 79 (March 2015) ; might be different version at the time of your analysis.
http://www.ensembl.org/info/data/ftp/index.html
then hit following commands:\n
gunzip *.gz
cat *.fa > "$sps"_v"$release"_dna.fa"
	exit 1;
fi

## download DNA fasta from latest version of ensemble;
# Note : danio_rerio link has diff than others sps

if [[ $sps == "dre" ]]; then
	echo -e "danio_rerio = "$sps"\n"
	echo -e "downloading your reference genome, please wait...\n"
	for chr in `seq 1 25` MT
	do
		if [[ $GRCh -eq 9 ]]; then
		$get $genome_dir/${chr}.$sps.tmp.fa.gz 2>/dev/null ftp://ftp.ensembl.org/pub/release-$release/fasta/danio_rerio/dna/Danio_rerio.Zv9.dna.chromosome.${chr}.fa.gz
			else
		$get $genome_dir/${chr}.$sps.tmp.fa.gz 2>/dev/null ftp://ftp.ensembl.org/pub/release-$release/fasta/danio_rerio/dna/Danio_rerio.GRCz$GRCh.dna.chromosome.${chr}.fa.gz
		fi
	done
	gunzip $genome_dir/*.gz
	cat $genome_dir/*.fa > $genome_dir/$sps"_v"$release.fa
	cat $genome_dir/$sps"_v"$release.fa | awk '{print $1}' | sed -e 's/^>/&chr/g' > $genome_dir/$sps.genome_v$release.fa
	rm -rf $genome_dir/*.tmp.fa $genome_dir/$sps"_v"$release.fa
	echo -e "Your Reference Genome FASTA has downloaded and formatted as prefer by the CBS-miRSeq pipeline. \n"

elif [[ $sps == "mmu" ]]; then
	echo -e "mus_musculus = "$sps"\n"
	echo -e "downloading your reference genome, please wait...\n"
	for chr in `seq 1 19` MT X Y
	do
		$get $genome_dir/${chr}.$sps.tmp.fa.gz 2>/dev/null  ftp://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/dna/Mus_musculus.GRCm"$GRCh".dna.chromosome.${chr}.fa.gz

	done
	gunzip $genome_dir/*.gz
	cat $genome_dir/*.fa > $genome_dir/$sps"_v"$release.fa
	cat $genome_dir/$sps"_v"$release.fa | awk '{print $1}' | sed -e 's/^>/&chr/g' > $genome_dir/$sps.genome_v$release.fa
	rm -rf $genome_dir/*.tmp.fa $genome_dir/$sps"_v"$release.fa
	echo -e "Your Reference Genome FASTA has downloaded and formatted as prefer by the CBS-miRSeq pipeline. \n"

elif [[ $sps == "rno" ]]; then
	echo -e "rattus_norvegicus = "$sps" "i.e. RAT"\n"
	echo -e "downloading your reference genome, please wait...\n"
	for chr in `seq 1 20` MT X Y
	do
		$get $genome_dir/${chr}.$sps.tmp.fa.gz 2>/dev/null  ftp://ftp.ensembl.org/pub/release-$release/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_"$GRCh".dna.chromosome.${chr}.fa.gz

	done
	gunzip $genome_dir/*.gz
	cat $genome_dir/*.fa > $genome_dir/$sps"_v"$release.fa
	cat $genome_dir/$sps"_v"$release.fa | awk '{print $1}' | sed -e 's/^>/&chr/g' > $genome_dir/$sps.genome_v$release.fa
	rm -rf $genome_dir/*.tmp.fa $genome_dir/$sps"_v"$release.fa
	echo -e "Your Reference Genome FASTA has downloaded and formatted as prefer by the CBS-miRSeq pipeline. \n"

elif [[ $sps == "hsa" ]]; then
	echo -e "human = "$sps"\n"
	echo -e "downloading your reference genome, please wait...\n"
	for chr in `seq 1 22` MT X Y
	do
		$get $genome_dir/${chr}.$sps.tmp.fa.gz 2>/dev/null ftp://ftp.ensembl.org/pub/release-$release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh"$GRCh".dna.chromosome.${chr}.fa.gz

	done
	gunzip $genome_dir/*.gz
	cat $genome_dir/*.fa > $genome_dir/$sps"_v"$release.fa
	cat $genome_dir/$sps"_v"$release.fa | awk '{print $1}' | sed -e 's/^>/&chr/g' > $genome_dir/$sps.genome_v$release.fa
	rm -rf $genome_dir/*.tmp.fa $genome_dir/$sps"_v"$release.fa
	echo -e "Your Reference Genome FASTA has downloaded and formatted as prefer by the CBS-miRSeq pipeline. \n"
else
	echo -e "SORRY: currently given "$sps" species not supported!! download it manually
and then run this script again to build bowtie index for CBS-miRSeq analysis !!\n"
exit 1;
fi


### checking if downloaded genome fasta exist
gfa=($genome_dir/*.genome*.fa)
if [[ ! -e "$gfa" ]]; then
	echo -e "ERROR !! Genome fasta has not downloaded, cant built inded."
	exit 1;
fi

#### building Index of downloaded genome
if [[ $index_color == "yes" ]]; then
	echo -e "Attention:\tbuilding COLOR-INDEX for your csfasta reads..please wait !!\n"
	bowtie-build -C $genome_dir/$sps.genome_v$release.fa $genome_dir/$sps.genome_v"$release"C > /dev/null
elif [ $index_color == "no" ]; then
	echo -e "Attention:\tbuilding BASE-SPACE index for your fastq reads..please wait !!\n"
	bowtie-build $genome_dir/$sps.genome_v$release.fa $genome_dir/$sps.genome_v"$release"B > /dev/null
else
	echo -e " index not build, please check inputs !!\n"
	exit 1;
fi

## create a index dir
mkdir -p $genome_dir/Index

## moving genome and index into created index dir
mv $genome_dir/*.{ebwt,fa,fasta} $genome_dir/Index 2>/dev/null
cp $genome_dir/Index/*.{fa,fasta} $genome_dir 2>/dev/null

##check if index has build successfully
bwt=($genome_dir/Index/*.ebwt)
if [[ ! -e "${bwt[0]}" ]]; then
	echo "ERROR : indexes are empty \n"
else
	echo -e "Index successfully built, now run CBS-miRseq module 1b for the further analysis..\n"
fi

echo -e "\nAnalysis finished:" `date`

## check if log exist
if [ -s "$logs/module1a.log" ]; then
	rm -f module1a.log
fi

) 2>&1) | tee -a module1a.log

mv -f module1a.log logs

echo -e "\n#For more info:log file also has created into the directory of CBS-miRSeq: logs/module1a.log\n"

exit $?

