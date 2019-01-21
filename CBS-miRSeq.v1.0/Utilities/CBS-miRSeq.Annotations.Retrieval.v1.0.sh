#!/bin/bash


## Date of creation: 18/04/2016
## This Module is a part of CBS-miRSeq help to fetch annotation based on release of
## Ensembl, RNAcentral and miRBase databases.
## Author:  Kesharwani RK; bioinforupesh2009.au@gmail.com

### check color package; before to proceeed.
if ! hash tput 2>/dev/null; then
	echo -e "tput does not exist; script aborting !!"
	echo -e "First install tput to run this script !!"
	echo -e "Hint: tput is part of the ncurses package"
	exit 0;
fi

echo -e '\0033\0143'
print_softInfo () {
echo ''
echo -e "#^^^^^^^^^^^^^^^Welcome to Annotation Retrieval Module^^^^^^^^^^^^^^^^^^^^^^^#"
echo -e "#		   Aim:	  To Retrieve Annotations		     	     #"
echo -e "#	 Requested Citation: Kesharwani RK et al.(2019)	     #"
echo -e "#		     Author:  Kesharwani RK; bioinforupesh2009.au@gmail.com	     #"
echo -e "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#"
echo ''
}
print_softInfo

echo -e "\n"
echo -e "\t\t~ ~ ~ ~ ~ CBS-miRSeq.Annotations.Retrieval.v1.0.sh~ ~ ~ ~ ~"
echo -e "\t\t\t~ ~ ~ ~ ~Annotations Retrieval~ ~ ~ ~ ~"
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

#Inputs: command line arguments
Output_dir="$1" ## OutPut directory
sps="$2" # three letter code of reference species
Ensembl_release="$3" ### ensembl release (two digit integer)
EnsemblAssemblyGRCh="$4" ## chromosomal assembly as numeric value (such as 38 or 36.67)
RNAcentral_release="$5" ## RNAcentral database release (two digit integer-> 1.0, 2.0, 3.0, 4.0, 5.0 ....)
RNAcentralAssemblyGRCh="$6" ## chromosomal assembly of RNAcentral as numeric value (such as 38 or 36.67)
miRBase_release="$7" ## miRBase release (two digit integer)
closely_related_sps="$8" ## closely related three letter code of reference species

## Probable closely related sps of reference species
## hsa==mmu
## mmu==hsa/rat
## dre==tni
## rno==mmu


## for checking length of sps
ref_cnt=$(echo -n "$sps" | wc -m)

print_USAGE()
{
echo -e "USAGE: bash ./CBS-miRSeq.module1a_v1.0.sh <Output_dir> <reference_species_3_letter_code> <Ensembl_release(two digit integer)> <EnsemblAssemblyGRCh> <RNAcentral_release(two digit integer)> <RNAcentralAssemblyGRCh> <miRBase_release(two digit integer)> <closely_related_sps(three_letter_code)>\n"
my_color_text "The example of Human Annotations:" green
my_color_text "bash ./CBS-miRSeq.module1a_v1.0.sh /.../annotation_dir/ hsa 78 38.78 2.0 38 21 mmu" magenta
my_color_text "The example of Mouse Annotations:" green
my_color_text "bash ./CBS-miRSeq.module1a_v1.0.sh /.../annotation_dir/ mmu 69 38.69 2.0 38 21 hsa" magenta
my_color_text "The example of Zebrafish Annotations:" green
echo "Zebrafish Genome Assembly 9"
my_color_text "bash ./CBS-miRSeq.module1a_v1.0.sh /.../annotation_dir/ dre 79 9 2.0 9 21 tni" magenta
echo "Zebrafish Genome Assembly 10"
my_color_text "bash ./CBS-miRSeq.module1a_v1.0.sh /.../annotation_dir/ dre 86 10 5.0 10 21 tni" magenta
my_color_text "The example of Rat Annotations:" green
echo "Rat Genome Assembly 5.0"
my_color_text "bash ./CBS-miRSeq.module1a_v1.0.sh /.../annotation_dir/ rno 79 5.0.79 2.0 5.0 21 mmu" magenta
echo "Rat Genome Assembly 6.0"
my_color_text "bash ./CBS-miRSeq.module1a_v1.0.sh /.../annotation_dir/ rno 86 6.0.86 5.0 6.0 21 mmu" magenta
echo -e "\n"
}

# checking input arguments
if [[ $# -ne 8 ]]; then
	echo -e "\n"
	my_color_text "#please supply all inputs.." cyan
	echo -e "\n"
	print_USAGE
	echo -e "#==============================================================================================================#"
	my_color_text "#NOTE: Probable Relation of version release of miRBase, Ensembl and RNAcentral database:" black
	my_color_text "mmu==miRbase v21==GenBank Assembly:GCA_000001635.2(mmu10)==Ensembl release v69(38.69)==RNAcentral 1.0/2.0" cyan
	my_color_text "hsa==miRbase v21==GenBank Assembly:GCA_000001405.15(hg38)= Ensembl release v78(38)==RNAcentral 1.0/2.0" cyan
	my_color_text "dre==miRbase v21==GenBank Assembly:GCA_000002035.2(Zv9)==Ensembl release v79(9)==RNAcentral 1.0/2.0" cyan
	my_color_text "rno==miRbase v21==GenBank Assembly:GCA_000001895.3(Rnor_5.0)==Ensembl release v79(5.0)==RNAcentral 1.0/2.0" cyan
	my_color_text "WARNNING: THE CBS-miRSeq pipeline IS TAKING NO GUARANTEE FOR THE RIGHT MAPPING FROM ONE DATABASE TO ANOTHER." red
	my_color_text "#Be careful with the input of Ensembl, RNAcentral and mRBase release and the chromosomal assembly" green
	my_color_text "visit for info: http://www.ensembl.org/info/data/ftp/index.html" magenta
	my_color_text "visit for info: ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/" magenta
	my_color_text "visit for info: http://www.mirbase.org/" magenta
	my_color_text "#NOTE: Probable closely related sps of reference species:" black
	my_color_text "mmu => hsa/rat" black
	my_color_text "hsa => mmu" black
	my_color_text "dre => tni" black
	my_color_text "rno => mmu/hsa" black
	my_color_text "#please supply all inputs.." cyan
	echo -e "\n"
	print_USAGE
	echo -e "#==============================================================================================================#"
	echo -e "\n"
	exit 0;
elif [[ ! -d "$Output_dir" ]]; then
    echo -e "\n#Input directory "$Output_dir" \ndoes not exist !!\n"
	my_color_text "#please supply all and correct inputs" cyan
	echo -e "\n"
	print_USAGE
	exit 1;
elif [[ $ref_cnt -ne 3 ]]; then
	echo -e "\n#Please provide 3 letter code of your reference species !!\n"
	echo -e "\nError(2)... \nplease supply all inputs.. \n"
	print_USAGE
	exit 1;
fi

## check perl script
extract_miR="./extract_miRNAs.pl"
if [[ ! -f $extract_miR ]]; then
	echo -e "Please insure that you are in current directory of CBS-miRSeq/Utilities.\n"
	echo -e "Sorry Job terminated.\n"
	exit 1;
fi

# # validation of Ensembl version release 
rel=$(echo -n "$Ensembl_release" | wc -c) ## counts numeric character which should be 2
if [[ "$rel" -ne 2 ]]; then
	echo -e "\n#ERROR: provided "$Ensembl_release" is not a correct integers,
please provide 2 digit integers only, make sure ensembl has this "$Ensembl_release""
	echo -e "\nError(4)... \nplease supply all inputs.. \n"
	print_USAGE
	exit 1;
else
	echo -e "\n"
	echo -e "#Analysis Species: "$sps""
	echo -e "#Closely Related sps: "$closely_related_sps""
	echo -e "#Ensembl Release: "$Ensembl_release""
	echo -e "#Ensembl Chromosomal Assembly: "$EnsemblAssemblyGRCh""
	echo -e "#RNAcentral Release: "$RNAcentral_release""
	echo -e "#RNAcentral Chromosomal Assembly: "$RNAcentralAssemblyGRCh""
	echo -e "#miRBase Release: "$miRBase_release""
	echo -e "#Output directory: "$Output_dir""
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

# checking supported sps for downloading corresponding genome from ensembl genome browser
if [[ $sps == "hsa" ]] || [[ $sps == "dre" ]] || [[ $sps == "mmu" ]] || [[ $sps == "rno" ]]; then
	echo -e " "
else
	echo -e "Currently supported reference Species to download..\n"
	echo -e "hsa	Human
dre	Zebrafish
mmu	Mouse
rno	Rat\n"
    echo -e "Given <"$sps"> species does not supported to download,
please download it manually and RUN this script again to build bowtie index.!!
Source to download genome from Ensembl(only Chr and MT; mandatory and preferred for CBS-miRSeq pipeline):
#Ensembl Release 79 (March 2015) ; might be different version at the time of your analysis.
http://www.ensembl.org/info/data/ftp/index.html
then hit following commands:\n
gunzip *.gz
cat *.fa > "$sps"_v"$Ensembl_release"_dna.fa"
	exit 1;
fi

## some warnings before to start
echo -e "\n"
echo -e "\t\t#### WARNING !! Before to Run this script ####"
echo -e "\n"
my_color_text "#Module is ready to process..." cyan
my_color_text "#Make sure internet connection works properly." red
my_color_text "#Also make certain that You have provided the correct arguments." green
echo -n "Continue ? (y/n) : "
read ans
echo -e "job running...\n"
echo -e "\n"
if [[ "${ans}" != "y" ]] && [[ "${ans}" != "Y" ]]; then
	clear
	echo -e "\n"
	my_color_text "Please note that any missing inputs/arguments can cause to fail annotation retrieval" red
	my_color_text "Thank you for using CBS-miRSeq pipeline." black
	echo -e "\n"
	exit 0;
fi

# # ################################ ANALYSIS BEGINS ######################################

## cleaning all existed files
rm -rf $Output_dir/*.{gz,fa,gff3,bed,gtf} 2>/dev/null
rm -rf $Output_dir/Annotations 2>/dev/null

### Retrieve Annotations From Ensembl database (gtf)

if [[ $sps == "hsa" ]]; then
	echo -e "human = "$sps"\n"
	$get $Output_dir/Homo_sapiens.GRCh$EnsemblAssemblyGRCh.Ensembl.gtf.gz 2>/dev/null ftp://ftp.ensembl.org/pub/release-$Ensembl_release/gtf/homo_sapiens/Homo_sapiens.GRCh$EnsemblAssemblyGRCh.gtf.gz
elif [[ $sps == "mmu" ]]; then
	echo -e "mus_musculus = "$sps"\n"
	$get $Output_dir/Mus_musculus.GRCm$EnsemblAssemblyGRCh.Ensembl.gtf.gz 2>/dev/null ftp://ftp.ensembl.org/pub/release-$Ensembl_release/gtf/mus_musculus/Mus_musculus.GRCm$EnsemblAssemblyGRCh.gtf.gz
elif [[ $sps == "rno" ]]; then
	echo -e "rattus_norvegicus = "$sps""/RAT"\n"
	$get $Output_dir/Rattus_norvegicus.Rnor_$EnsemblAssemblyGRCh.Ensembl.gtf.gz 2>/dev/null ftp://ftp.ensembl.org/pub/release-$Ensembl_release/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_$EnsemblAssemblyGRCh.gtf.gz
fi

## zebrafish release check

if [[ $sps == "dre" ]]; then
	if [[ $EnsemblAssemblyGRCh -eq "10" ]]; then
		echo -e "danio_rerio = "$sps"\n"
		$get $Output_dir/Danio_rerio.GRCz$EnsemblAssemblyGRCh.Ensembl.gtf.gz 2>/dev/null ftp://ftp.ensembl.org/pub/release-$Ensembl_release/gtf/danio_rerio/Danio_rerio.GRCz$EnsemblAssemblyGRCh.$Ensembl_release.gtf.gz
	elif [[ $EnsemblAssemblyGRCh -eq "9" ]]; then
		echo -e "danio_rerio = "$sps"\n"
		$get $Output_dir/Danio_rerio.GRCz$EnsemblAssemblyGRCh.Ensembl.gtf.gz 2>/dev/null ftp://ftp.ensembl.org/pub/release-$Ensembl_release/gtf/danio_rerio/Danio_rerio.Zv$EnsemblAssemblyGRCh.$Ensembl_release.gtf.gz
	else
		echo -e "Program does not know Ensembl_release for the Zerafish,
		Please download it manually."
	fi
fi
# ## check files
if [ -s $Output_dir/*Ensembl.gtf.gz ]; then
	echo -e "Retrieval of Ensembl annotation is done.\n"
else
	my_color_text "ERROR !! Ensembl gtf is empty." red
	echo "Please recheck input of your Ensembl Assembly and version release."
	my_color_text "visit for info: http://www.ensembl.org/info/data/ftp/index.html" magenta
	echo "OR download it manually."
	echo -e "\n"
	print_USAGE
	exit 1;

fi

### Retrieve Annotation From RNAcentral database (bed)

if [[ $sps == "hsa" ]]; then
	$get $Output_dir/Homo_sapiens.GRCh$RNAcentralAssemblyGRCh.RNAcentral.bed.gz 2>/dev/null ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/$RNAcentral_release/genome_coordinates/Homo_sapiens.GRCh$RNAcentralAssemblyGRCh.bed.gz
elif [[ $sps == "mmu" ]]; then
	$get $Output_dir/Mus_musculus.GRCm$RNAcentralAssemblyGRCh.RNAcentral.bed.gz 2>/dev/null ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/$RNAcentral_release/genome_coordinates/Mus_musculus.GRCm$RNAcentralAssemblyGRCh.bed.gz
elif [[ $sps == "rno" ]]; then
	$get $Output_dir/Rattus_norvegicus.Rnor_$RNAcentralAssemblyGRCh.RNAcentral.bed.gz 2>/dev/null ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/$RNAcentral_release/genome_coordinates/Rattus_norvegicus.Rnor_$RNAcentralAssemblyGRCh.bed.gz
fi

## zebrafish release check
if [[ $sps == "dre" ]]; then
	if [[ $EnsemblAssemblyGRCh -eq "10" ]]; then
		echo -e "danio_rerio = "$sps"\n"
		$get $Output_dir/Danio_rerio.GRCz$RNAcentralAssemblyGRCh.RNAcentral.bed.gz 2>/dev/null ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/$RNAcentral_release/genome_coordinates/Danio_rerio.GRCz$RNAcentralAssemblyGRCh.bed.gz
	elif [[ $EnsemblAssemblyGRCh -eq "9" ]]; then
		echo -e "danio_rerio = "$sps"\n"
		$get $Output_dir/Danio_rerio.GRCz$RNAcentralAssemblyGRCh.RNAcentral.bed.gz 2>/dev/null ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/$RNAcentral_release/genome_coordinates/Danio_rerio.Zv$RNAcentralAssemblyGRCh.bed.gz
	else
		echo -e "Program does not know RNAcentral Assembly for the Zerafish,
		Please download it manually."
	fi
fi

# ## check files
if [ -s $Output_dir/*RNAcentral.bed.gz ]; then
	echo -e "Retrieval of RNAcentral annotation is done.\n"
else
	my_color_text "ERROR !! RNAcentral bed is empty." red
	echo "Please recheck input of your RNAcentral Assembly and version release."
	my_color_text "visit for info: ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/" magenta
	echo -e "\n"
	echo "OR download it manually."
	print_USAGE
	exit 1;
fi

# ### Retrieve Annotations From miRbase (gff3, mature fasta, closely related sps mature fasta and precuros fasta)

$get $Output_dir/$sps.miRBase_v$miRBase_release.gff3 2>/dev/null ftp://mirbase.org/pub/mirbase/$miRBase_release/genomes/$sps.gff3
# ## check files
if [ ! -s $Output_dir/*miRBase_v$miRBase_release.gff3 ]; then
	my_color_text "ERROR !! miRBase gff3 is empty." red
	echo "Please recheck input of your miRBase version release."
	my_color_text "visit for info: http://www.mirbase.org/" magenta
	echo -e "\n"
	exit 1;
fi

$get $Output_dir/All.sps.mature_v$miRBase_release.fa.gz 2>/dev/null ftp://mirbase.org/pub/mirbase/$miRBase_release/mature.fa.gz
# ## check files
if [ ! -s $Output_dir/All.sps.mature_v$miRBase_release.fa.gz ]; then
	echo "ERROR !! miRBase mature.fa is empty."
	echo "Please recheck input of your miRBase version release."
	my_color_text "visit for info: http://www.mirbase.org/" magenta
	echo -e "\n"
	exit 1;
fi

$get $Output_dir/All.sps.hairpin_v$miRBase_release.fa.gz 2>/dev/null ftp://mirbase.org/pub/mirbase/$miRBase_release/hairpin.fa.gz
# ## check files
if [ ! -s $Output_dir/All.sps.hairpin_v$miRBase_release.fa.gz ]; then
	my_color_text "ERROR !! miRBase hairpin.fa is empty." red
	echo "Please recheck input of your miRBase version release."
	my_color_text "visit for info: http://www.mirbase.org/" magenta
	echo -e "\n"
	exit 1;
fi

## Uncompressed all files
gunzip $Output_dir/*.gz

### Retrieve Annotation From miRBase database (fetch mature, precursor and Related sps fasta)
perl $extract_miR $Output_dir/All.sps.mature_v$miRBase_release.fa $sps > $Output_dir/$sps.mature_v$miRBase_release.fa
perl $extract_miR $Output_dir/All.sps.hairpin_v$miRBase_release.fa $sps > $Output_dir/$sps.precursor_v$miRBase_release.fa
## fetching closely related sps fasta
perl $extract_miR $Output_dir/All.sps.mature_v$miRBase_release.fa $closely_related_sps > $Output_dir/$closely_related_sps.Related.mature_v$miRBase_release.fa
echo -e "Retrieval of miRBase annotation is done.\n"

## put all annotations into dir
mkdir -p $Output_dir/Annotations
mv -f $Output_dir/*.{fa,gff3,bed,gtf} $Output_dir/Annotations 2>/dev/null

my_color_text "All annotation successfully fetched." cyan
my_color_text "Now you may proceed with your analyses using modules of CBS-miRSeq." green
my_color_text "Thank you for using the CBS-miRSeq pipeline." green
echo -e "\n"

