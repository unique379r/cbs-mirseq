#!/bin/bash

clear

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
# # USAGE: Rscript bash ./Install.CBS-miRSeq.dependencies.v1.0.sh
# # Author: Kesharwani RK; Email: bioinforupesh2009.au@gmail.com
# # version: 1.0
# # Copyright (c) 2019 Kesharwani RK

echo -e "#============================= Dependencies installation ================================================#"
echo -e "## Description: This script aims to install dependencies software and tools for the CBS-miRSeq pipeline"
echo -e "## Copyright (c) 2019 Kesharwani RK"
echo -e "## version: 1.0"
echo -e "## Author: Kesharwani RK; Email: bioinforupesh2009.au@gmail.com"
echo -e "#========================================================================================================#"

### check color package; before to proceeed.
if ! hash tput >/dev/null 2>&1; then
	echo -e "tput does not exist; script aborting !!"
	echo -e "First install tput to run this script !!"
	exit 0;
fi

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

## exit program with msg
set -e
die() {
    echo -e "FATAL ERROR: $* (status $?)" 1>&2
    exit 1;
}

## some warnings before to start
echo -e "\n"
echo -e "\t\t\t\t#### Before to Run this script ####"
echo -e "\n"
my_color_text "#Make sure you have successfully installed miRDeep2 (v2.0.0.5) along with their dependencies." blue
my_color_text "#Ignore bowtie installation if you already have installed during miRDeep2 installtion." blue
my_color_text "#Make sure internet connection works properly under super-user or root privileges." cyan
my_color_text "#This script must be run with super-user (use su -) permission." red
my_color_text "#Example to enter super user as follows:" blue
echo -e "$ su -"
echo -e "Password: <TYPE ROOT PASSWORD>"
echo -e "# bash ./install.CBS-miRSeq.dependencies.v1.0.sh"
my_color_text "#Any dependency already installed on the system has to be available from the \$PATH (Global) environment variable." magenta
echo -n "Continue ? (y/n) : "
read ans
if [[ "${ans}" != "y" ]] && [[ "${ans}" != "Y" ]]; then
	echo -e "\n"
	clear
	my_color_text "Please note that without dependencies, CBS-miRSeq is not able to proceed analysis !!" red
	my_color_text "Thank you for using CBS-miRSeq pipeline." blue
	exit 0;
fi

# # recheck sudo access
if [[ $EUID -ne 0 ]]; then
	my_color_text "ERROR !! Sorry, you are not root." red
	echo -e "Please contact your system administrator.\n"
	echo -e "\n"
	exit 1;
fi

echo -e "\n"
my_color_text "#Checking dependencies ..." magenta
echo -e "\n"

# #check OS (Unix/Linux or Mac)
os=`uname`;

# # get the right download program based on OS
if [[ "$os" = "Darwin" ]]; then
	# use curl as the download program
	get="curl -L -o"
else
	# use wget as the download program
	get="wget --no-check-certificate -O"
fi

########### Checking all dependencies software if installed; if no install it ###################

##first thing first to check miRDeep2
### check miRDeep2 package; before proceeed.
## To check if program is executable and set in PATH
##if [[ ! -x "$(command -v miRDeep2.pl)" ]]; then
	##echo -e "miRDeep2 is not installed, please install first prior to use this installation.\n" >&2
##else
	##echo -e "miRDeep installed in your system.\n"
##fi

###############################
#### system utilities check ###
###############################

### check unzip package; before proceeed.
if ! hash unzip >/dev/null 2>&1; then
	my_color_text "Error !! Can not proceed without unzip, please install:" red
	my_color_text "zlib1g-dev" cyan
	echo -e "\n"
	exit 1;
fi

### check gcc package; before proceeed.
if ! hash gcc >/dev/null 2>&1; then
	my_color_text "Error !! Can not proceed without GCC, please install:" red
	my_color_text "gcc" cyan
	echo -e "\n"
	exit 1;
fi

## C++ libraries
if [[ ! -e /usr/include/zlib.h ]]; then
    my_color_text "Error !! zlib.h library is missing. please install:" red
	my_color_text "zlib1g-dev" cyan
	echo -e "\n"
    exit 1;
fi

## check libncurses5-dev
if [[ ! -e /usr/include/ncurses.h ]]; then
    my_color_text "Error !! ncurses.h not found. please install:" red
	my_color_text "libncurses5-dev" cyan
	echo -e "\n"
    exit 1;
fi

### check g++ package; before proceeed.
if ! hash g++ >/dev/null 2>&1; then
	my_color_text "Error !! Can not proceed without g++, please install:" red
	my_color_text "http://gcc.gnu.org/" cyan
	echo -e "\n"
	exit 1;
fi


### check make package; before proceeed.
if ! hash make >/dev/null 2>&1; then
	my_color_text "Error !! Can not proceed without make, please install:" red
	my_color_text "ftp://ftp.gnu.org/gnu/make/" cyan
	echo -e "\n"
	exit 1;
fi

### check make package; before proceeed.
if ! hash automake >/dev/null 2>&1; then
	my_color_text "Error !! Can not proceed without automake, please install:" red
	my_color_text "http://www.gnu.org/software/automake/#downloading" cyan
	echo -e "\n"
	exit 1;
fi

### check awk package; before proceeed.
if ! hash awk >/dev/null 2>&1; then
    my_color_text "Can not proceed without awk, please install:" red
	my_color_text "gawk" cyan
	echo -e "\n"
	exit 1;
fi

### check perl package; before proceeed.
if ! hash perl >/dev/null 2>&1; then
    my_color_text "Can not proceed without perl, please install:" red
    my_color_text "perl" cyan
	echo -e "\n"
	exit 1;
fi

### check python package; before proceeed.
if ! hash python >/dev/null 2>&1; then
    my_color_text "Can not proceed without python, please install:" red
    my_color_text "python >= 2.7" cyan
	echo -e "\n"
	exit 1;
fi

### check java package; before proceeed.
if ! hash java >/dev/null 2>&1; then
    my_color_text "Can not proceed without java, please install:" red
    my_color_text "java v1.6-v1.8 JREs from Oracle" cyan
	echo -e "\n"
	exit 1;
fi


################ Set preferred installation directory ###################

echo "Where should missing software prerequisites be INSTALLED ? (Please give absolute path) "
read ans
#ans=${ans:-$PREFIX_BIN}
PREFIX_BIN=$ans
if [[ ! -d $PREFIX_BIN ]]; then
    echo "Directory $PREFIX_BIN does not exist!"
    echo -n "Do you want to create $PREFIX_BIN folder ? (y/n)  : "
    read reply
    if [[ "${reply}" = "y" || "${reply}" = "Y" ]]; then
	mkdir -p $PREFIX_BIN/bin
    else
	die "Must specify a directory to install required software!"
    fi
fi
if [[ ! -w $PREFIX_BIN ]]; then
    die "Cannot write to directory $PREFIX_BIN."
fi

echo -e "\n"
echo -e "\n"
################ Set preferred source directory ###################
echo "Where should missing software be DOWNLOAD ? (Please give absolute path) "
read ans
source=$ans
if [[ ! -d $source ]]; then
    echo "Directory $source does not exist!"
    echo -n "Do you want to create $source folder ? (y/n)  : "
    read reply
    if [[ "${reply}" = "y" || "${reply}" = "Y" ]]; then
	mkdir -p $source
    else
	die "Must specify a directory to download and install required software!"
    fi
fi

if [[ ! -w $source ]]; then
    die "Cannot write to directory $source."
fi


## recreate
mkdir -p $source


## create user defined PREFIX_BIN
mkdir -p $PREFIX_BIN/bin
PREFIX_BIN=$PREFIX_BIN/bin

#########################
#### Install Programs ###
#########################
clear
echo -e "Processing...."
##### Bowtie #####
## bowtie v1
# is already installed?
if hash bowtie >/dev/null 2>&1; then
	echo -e "\n"
	bowtie_version=`bowtie --version | head -1`;
	my_color_text "#bowtie ($bowtie_version) appear to have already installed !" blue
else
	echo -e "\n"
	my_color_text "#Bowtie Aligner appear to have NOT installed !" blue
fi
echo -n "#Would you like to install/re-install Bowtie ? (y/n)  : "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading bowtie, please wait...."
	$get $source/bowtie-1.0.0-linux-x86_64.zip 2>/dev/null http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip/
	#$get $source/bowtie-0.12.9-linux-x86_64.zip 2>/dev/null http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.9/bowtie-0.12.9-linux-x86_64.zip/download
	cd $source/
	unzip bowtie-1.0.0-linux-x86_64.zip
	#unzip bowtie-0.12.9-linux-x86_64.zip
	cd bowtie-1.0.0
	#cd bowtie-0.12.9
	echo -e "\ncopying excutable into Path.\n"
	cp -rav bowtie bowtie-build bowtie-inspect $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/bowtie ]]; then
	my_color_text "bowtie is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without bowtie" red
	echo -e "\n"
else
	my_color_text "Can not proceed without bowtie, please install it manually and re-launch this script" red
	exit 1;
fi

### BEDTools #####
# is already installed?
if hash bedtools >/dev/null 2>&1; then
	echo -e "\n"
	bedtools_version=`bedtools --version`;
	my_color_text "#BEDTools ($bedtools_version) appear to have already installed !" blue
else
	echo -e "\n"
	my_color_text "#BEDTools appear to have NOT installed !" blue
fi

echo -n "#Would you like to install/re-install BEDTools ? (y/n)	 : "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading bedtools, please wait...."
	$get $source/bedtools-2.23.0.tar.gz 2>/dev/null http://github.com/arq5x/bedtools2/releases/download/v2.23.0/bedtools-2.23.0.tar.gz
	cd $source
	tar -zxvf bedtools-2.23.0.tar.gz
	cd bedtools2
	echo -e "\ncopying excutable into Path.\n"
	make clean
	make all
	cp -rav bin/* $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/bedtools ]]; then
	my_color_text "bedtools is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without bedtools" red
	echo -e "\n"
else
	my_color_text "Can not proceed without bedtools, please install it manually and re-launch this script" red
	exit 1;
fi

# # ###### samtools ######
# is already installed?
if hash samtools >/dev/null 2>&1; then
	echo -e "\n"
	my_color_text "#samtools appear to have already installed !" blue
else
	echo -e "\n"
	my_color_text "#samtools appear to have NOT installed !" blue
fi
echo -n "#Would you like to install/re-install samtools ? (y/n)	 : "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading samtools, please wait...."
	$get $source/samtools-0.1.18.tar.bz2 2>/dev/null http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download
	#$get $source/samtools-0.1.18.tar.bz2 2>/dev/null http://netix.dl.sourceforge.net/project/samtools/samtools/0.1.8/samtools-0.1.8.tar.bz2
	cd $source
	tar -xvjpf samtools-0.1.18.tar.bz2
	cd samtools-0.1.18
	echo -e "\ncopying excutable into Path.\n"
	make
	cp -rav samtools $PREFIX_BIN
	cp -rav misc/md5fa misc/md5sum-lite misc/wgsim misc/*.pl bcftools/bcftools $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/samtools ]]; then
	my_color_text "samtools is installed successfully." cyan
	# # Add binary directory to path globally
	echo 'export PATH=$PATH:'$PREFIX_BIN >> ~/.bashrc
	echo 'PATH=$PATH:'$PREFIX_BIN >> ~/.bash_profile
	## reloading .bashrc
	source ~/.bashrc
	source ~/.bash_profile
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without samtools" red
	echo -e "\n"
else
	my_color_text "Can not proceed without samtools, please install it manually and re-launch this script" red
	exit 1;
fi

# # ###### samstat ######
# is already installed?
if hash samstat >/dev/null 2>&1; then
	echo -e "\n"
	samstat_version=`samstat --version`
	my_color_text "#samstat ($samstat_version) appear to have already installed !" blue
else
	echo -e "\n"
	my_color_text "#samstat appear to have NOT installed !" blue
fi

echo -n "#Would you like to install/re-install samstat ? (y/n)	: "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading samstat, please wait...."
	$get $source/samstat-1.5.1.tar.gz 2>/dev/null http://sourceforge.net/projects/samstat/files/samstat-1.5.1.tar.gz/download
	cd $source
	tar -zxvf samstat-1.5.1.tar.gz
	cd samstat-1.5.1
	echo -e "\ncopying excutable into Path.\n"
	./configure && make && make install
	cp -v src/samstat $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/samstat ]]; then
	my_color_text "samstat is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without samstat" red
	echo -e "\n"
else
	my_color_text "Can not proceed without samstat, please install it manually and re-launch this script" red
	exit 1;
fi


### FASTQC #####
# is already installed?
if hash fastqc >/dev/null 2>&1; then
	echo -e "\n"
	fastqc_version=`fastqc -v`
	my_color_text "#fastqc ($fastqc_version) appear to have already installed !" blue
else
	echo -e "\n"
	my_color_text "#fastqc appear to have NOT installed !" blue
fi

echo -n "Would you like to install/re-install fastqc? (y/n)  : "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading FastQC, please wait...."
	$get $source/fastqc_v0.10.1.zip 2>/dev/null http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip
	cd $source
	unzip fastqc_v0.10.1.zip
	cd FastQC
	echo -e "\ncopying excutable into Path.\n"
	chmod -R 755 fastqc
	cp -av fastqc $PREFIX_BIN
	#sudo ln -sf $source/FastQC/fastqc /usr/loca/bin/fastqc
	#sudo ln -sf $source/FastQC/fastqc $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/fastqc ]]; then
	my_color_text "fastqc is installed successfully." cyan
	sudo ln -sf $source/FastQC/fastqc /usr/local/bin/fastqc
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without fastqc" red
	echo -e "\n"
else
	my_color_text "Can not proceed without fastqc, please install it manually and re-launch this script" red
	exit 1;
fi

### FASTX #####
echo;
# is already installed?
if hash fastq_quality_filter >/dev/null 2>&1; then
	echo -e "\n"
	fastx_version=`fastq_quality_filter -h | grep -i Part | awk '{print $3,$4,$5}'`
	my_color_text "#fastx ($fastx_version) appear to have already installed !" blue
else
	echo -e "\n"
	my_color_text "#fastx appear to have NOT installed !" blue
fi
echo -n "Would you like to install/re-install fastx? (y/n)	: "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading FastX, please wait...."
	$get $source/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 2>/dev/null http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
	cd $source
	tar -xvjpf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
	#echo "Renaming bin to fastx_toolkit_0.0.13"
	mv ./bin ./fastx_toolkit_0.0.13
	#echo entering fastx_toolkit_0.0.13
	#cd fastx_toolkit_0.0.13
	#ls
	echo -e "\ncopying excutable into Path.\n"
	##./configure && make && make install
	cp -av fastx_toolkit_0.0.13/* $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/fastq_quality_filter ]]; then
	my_color_text "fastx is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without fastx" red
	echo -e "\n"
else
	my_color_text "Can not proceed without fastx, please install it manually and re-launch this script" red
	exit 1;
fi

### cutadapt #####
# is already installed?
if hash cutadapt >/dev/null 2>&1; then
	echo -e "\n"
	cutadapt_version=`cutadapt --version`
  my_color_text "#cutadapt v"$cutadapt_version" appear to have already installed !" blue
	else
	echo -e "\n"
	my_color_text "#cutadapt appear to have NOT installed !" blue
fi

echo -n "Would you like to install/re-install cutadapt? (y/n)  : "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading cutadapt, please wait...."
	$get $source/cutadapt-1.3.tar.gz 2>/dev/null https://www.dropbox.com/s/l9u22kr8zbgi3ky/cutadapt-1.3.tar.gz?dl=1
	##mkdir -p $source/cutadpat_latest
	##git clone --depth=1 https://github.com/marcelm/cutadapt.git $source/cutadpat_latest
	cd $source
	tar -zxvf cutadapt-1.3.tar.gz
	cd cutadapt-1.3
	echo -e "\ncopying excutable into Path.\n"
	python setup.py build
	python setup.py install
	cp -av bin/cutadapt $PREFIX_BIN
	cd ..
	
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/cutadapt ]]; then
	my_color_text "cutadapt is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without cutadapt" red
	echo -e "\n"
else
	my_color_text "Can not proceed without cutadapt, please install it manually and re-launch this script" red
	exit 1;
fi
### ## featureCounts v1.4.6 #####
# is already installed?
if hash featureCounts >/dev/null 2>&1; then
	echo -e "\n"
	featureCounts_version=`featureCounts -v`
	my_color_text "# featureCounts "$featureCounts_version" appear to have already installed !" blue
	else
	echo -e "\n"
	my_color_text "#featureCounts appear to have NOT installed !" blue
fi

echo -n "Would you like to install/re-install featureCounts? (y/n)	: "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading featureCounts, please wait...."
	$get $source/subread-1.4.6-p1-Linux-x86_64.tar.gz 2>/dev/null http://sourceforge.net/projects/subread/files/subread-1.4.6-p1/subread-1.4.6-p1-Linux-x86_64.tar.gz/download
	cd $source
	tar -zxvf subread-1.4.6-p1-Linux-x86_64.tar.gz
	cd subread-1.4.6-p1-Linux-x86_64 ## its binary file no need of installation
	echo -e "\ncopying excutable into Path.\n"
	cp -av bin/featureCounts $PREFIX_BIN
	cd ..
	cd ..
  # installed=0;
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/featureCounts ]]; then
	my_color_text "featureCounts is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without featureCounts" red
	echo -e "\n"
else
	my_color_text "Can not proceed without featureCounts, please install it manually and re-launch this script" red
	exit 1;
fi

### RNAhybrid v2.1.1 #####
# is already installed?
if hash RNAhybrid >/dev/null 2>&1; then
	echo -e "\n"
	my_color_text	"#RNAhybrid appear to have already installed !" blue
	else
	echo -e "\n"
	my_color_text "#RNAhybrid appear to have NOT installed !" blue
fi
echo -n "Would you like to install/re-install RNAhybrid? (y/n)	: "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading RNAhybrid, please wait...." ## problem to downloading from main site as he asked for choice of file types
	## however, here i uploaded file into dropbox to download thew file without any cause
	$get $source/RNAhybrid-2.1.1-src.tar.gz 2>/dev/null https://www.dropbox.com/s/qvpb4yu7agjitby/RNAhybrid-2.1.1-src.tar.gz?dl=1
	cd $source
	tar -xvf RNAhybrid-2.1.1-src.tar.gz
	cd RNAhybrid-2.1.1
	echo -e "\ncopying excutable into Path.\n"
	./configure && make && make install
	cp -av src/RNAcalibrate src/RNAhybrid src/RNAeffective $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/RNAhybrid ]]; then
	my_color_text "RNAhybrid is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without RNAhybrid" red
	echo -e "\n"
else
	my_color_text "Can not proceed without RNAhybrid, please install it manually and re-launch this script" red
	exit 1;
fi

### miRanda 3.3.a #####
# is already installed?
if hash miranda >/dev/null 2>&1; then
	echo -e "\n"
	miranda_version=`miranda -v | grep -i 'microRNA Target' | awk '{print $1,$2}'`
	my_color_text	"#miranda ($miranda_version) appear to have already installed !" blue
else
	echo -e "\n"
	my_color_text "#miranda appear to have NOT installed !" blue
fi

echo -n "Would you like to install/re-install miranda? (y/n)  : "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading miRanda, please wait...."
	$get $source/miRanda-aug2010.tar.gz 2>/dev/null http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz
	cd $source
	tar -zxvf miRanda-aug2010.tar.gz
	cd miRanda-3.3a
	echo -e "\ncopying excutable into Path.\n"
	./configure && make && make install
	cp -av src/miranda $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/miranda ]]; then
	my_color_text "miRanda is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without miranda" red
	echo -e "\n"
else
	my_color_text "Can not proceed without miRanda, please install it manually and re-launch this script" red
	exit 1;
fi

##### solid2fastq from bfast #####
# is already installed?
if hash solid2fastq >/dev/null 2>&1; then
	echo -e "\n"
	my_color_text "#solid2fastq of bfast appear to have already installed.!" blue
	else
	echo -e "\n"
	my_color_text "#solid2fastq of bfast appear to have NOT installed !" blue
fi
echo -n "#Would you like to install/re-install solid2fastq of bfast ? (y/n)  : "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -e "\nDownloading bfast, please wait...."
	$get $source/bfast-0.7.0a.tar.gz 2>/dev/null http://sourceforge.net/projects/bfast/files/bfast/0.7.0/bfast-0.6.5a.tar.gz/download
	cd $source
	tar -zxvf bfast-0.6.5a.tar.gz
	cd bfast-0.6.5a
	echo -e "\ncopying excutable into Path.\n"
	sh autogen.sh
	./configure
	make
	sudo make install
	cp -av scripts/solid2fastq $PREFIX_BIN
	cd ..
	cd ..
fi
clear
### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/solid2fastq ]]; then
	my_color_text "solid2fastq from bfast is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without solid2fastq of bfast" red
	echo -e "\n"
else
	my_color_text "Can not proceed without solid2fastq of bfast, please install it manually and re-launch this script" red
	exit 1;
fi


##### R/BioConductor #####
# is already installed?
if hash R >/dev/null 2>&1; then
	echo -e "\n"
	R_version=`R --version | head -1 | awk '{print $2,$3}'`
	my_color_text "#R $R_version appear to have already installed !" blue
else
	echo -e "\n"
	my_color_text "#R appear to have NOT installed !" blue
fi
echo -n "Would you like to install/re-install R? (y/n)	: "
read ans
if [[ "${ans}" = "y" || "${ans}" = "Y" ]]; then
	echo -n "#Please provide a version release (three digit like 3.3.2) that need to be installed	: "
	read var
	echo -e "\nDownloading R, please wait...."
    $get $source/R-${var}.tar.gz 2>/dev/null https://cran.r-project.org/src/base/R-3/R-${var}.tar.gz
	## to check if downloading failed
	if [[ ! -s $source/R-${var}.tar.gz ]] ; then
		rm -rf $source/R-${var}.tar.gz
		echo -e "#Please confirm the version release of R."
		echo -e "#Visit for more info: https://cran.r-project.org/src/base/R-3/"
		echo -e "OR, Install it manually and re-launch this script."
		exit 1;
	fi
	cd $source
	tar -zxvf R-${var}.tar.gz
	cd R-${var}
	./configure --with-x=no && make
	echo -e "\ncopying excutable into Path.\n"
	cp -av bin/R bin/Rscript $PREFIX_BIN
	cd ..
	cd ..
fi
clear

## If error found then probably need 'flavours' of the libcurl4 development package to install R (≥ 3.0.0)
##libcurl4-gnutls-dev
##libcurl4-nss-dev
##libcurl4-openssl-dev (prefered)
## https://cran.r-project.org/web/packages/curl/index.html
## sudo apt-get install libcurl4-openssl-dev


### double check if excutable file has copied..
if [[ -x $PREFIX_BIN/Rscript ]]; then
	my_color_text "R is installed successfully." cyan
	echo -e "\n"
elif [[ "${ans}" = "n" || "${ans}" = "N" || "${ans}" = "" ]]; then
	my_color_text "Please note that CBS-miRSeq can not proceed without R" red
	echo -e "\n"
else
	my_color_text "Can not proceed without R, please install it manually and re-launch this script" red
	exit 1;
fi

# # ##Apply the permission to all the files under a directory recursively
##chmod -R 777 "$PREFIX_BIN/"*
echo 'PATH=$PATH:'$PREFIX_BIN >> ~/.bash_profile
echo 'export PATH=$PATH:'$PREFIX_BIN >> ~/.bashrc
## reload terminal
source ~/.bashrc
source ~/.bash_profile
clear
################ End of the Installation ###################
echo -e "\n"
my_color_text "Installation done !" magenta
my_color_text "CBS-miRSeq dependencies are installed correctly, now you may proceed with your analysis." yellow
my_color_text "Thank you for using CBS-miRSeq pipeline." green
echo -e "\n"

# # ##list Dependencies
cbs_dep="bowtie bedtools samtools samstat solid2fastq fastqc(with_option-h) fastq_quality_filter(with_option-h) cutadapt(with_option-h) featureCounts RNAhybrid(with_option-h) miranda(with_option-h) RNAfold(with_option-h) randfold make_html.pl miRDeep2.pl R(with_option--version) "
length=$(echo $cbs_dep | wc -w)
my_color_text "Recommended step !" blue
echo -e "To varify if everything is installed properly type:\n"
for dep in $cbs_dep; do
    echo "${dep}"
done
echo -e "\n"

echo -e "#PLEASE INGNORE THE WARNING (Unescaped left brace in regex is deprecated, passed through in regex; marked by) MAY FOUND BY make_html.pl"

my_color_text "Please close and reopen your terminal to run above listed commands.. !" blue
my_color_text "Also Please Run R packages installation script in order to install required R packages!" blue
echo -e "\n"
echo -e "USAGE:"
echo -e "\n"
echo -e "# Rscript ./CBS-miRSeq.Required.Packages.R"
echo -e "# If required, Please make change of permission/ownership of:" $source 
echo -e "# sudo chown -R UserName:UserName"$source
echo -e "\n"
echo -e "Done.\n"

###
exec bash --login



