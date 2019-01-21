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
# # Author: Kesharwani RK email: bioinforupesh2009.au@gmail.com
# # version: 1.0
# # Copyright (c) 2019 Kesharwani RK

echo -e "#================ Welcome to System Dependencies installation Script =====================================#"
echo -e "## Description: This script aims to install System Dependencies and lacking packages"
echo -e "## Copyright (c) 2019 Kesharwani RK"
echo -e "## version: 1.0"
echo -e "## Author: Kesharwani RK; Email: bioinforupesh2009.au@gmail.com"
echo -e "#=========================================================================================================#"

### USAGE: bash ./CBS-miRSeq-SystemPackagesInstall.v1.0.sh


### check color package; before proceeed.
if ! hash tput 2>/dev/null; then
	echo -e "tput does not exist; script aborting !!"
	echo -e "First install tput package to run this script !!"
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

## some warnings before to start
echo -e "\n"
echo -e "\t\t\t\t#### Before to Run this script ####"
echo -e "\n"
my_color_text "#Make sure internet connection works properly under super-user or root privileges." cyan
my_color_text "#This script must be run with super-user (use su -) permission." red
my_color_text "#Example to enter super user as follows:" blue
echo -e "$ su -"
echo -e "Password: <TYPE ROOT PASSWORD>"
echo -e "# bash ./CBS-miRSeq-SystemPackagesInstall.v1.0.sh"
echo -n "Continue ? (y/n) : "
read ans
if [[ "${ans}" != "y" ]] && [[ "${ans}" != "Y" ]]; then
	echo -e "\n"
	clear
	my_color_text "Please note that without system packages, dependencies script CBS-miRSeq may not be able to run successfully !!" red
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
my_color_text "#Checking missing dependencies on your system ..." magenta
echo -e "\n"


## define DEPENDENCIES into an array
declare -a DEPENDENCIES=("unzip" "gzip" "wget" "gcc" "bzip2" "g++" "git-all" "make" "automake" "perl" "python" "java")
# What dependencies are missing?
PKGSTOINSTALL=""
for i in "${DEPENDENCIES[@]}"; do
	# # for Ubuntu like system
	if hash dpkg &> /dev/null; then
		if [[ `dpkg -l $i` == *'is not installed'* ]]; then
			PKGSTOINSTALL=$PKGSTOINSTALL" "${DEPENDENCIES[$i]}
		fi
	# # # # # OpenSuse, Mandriva, Fedora, CentOs, ecc. (with rpm)
	elif hash rpm &> /dev/null; then
		if [[ `rpm -q $i` == *'is not installed'* ]]; then
			PKGSTOINSTALL=$PKGSTOINSTALL" "$i
		fi
	# ArchLinux (with pacman)
	elif hash pacman &> /dev/null; then
		if [[ `pacman -Qqe $i` == *'is not installed'* ]]; then
			PKGSTOINSTALL=$PKGSTOINSTALL" "${DEPENDENCIES[$i]}
		fi
	# If its impossible to determine if there are missing dependencies, mark all as missing
	else
		echo "last"
		PKGSTOINSTALL=$PKGSTOINSTALL" "$i
		echo "finish"
	fi
done

# # echo -e "These DEPENDENCIES need to be installed:\n"
# # echo -e $PKGSTOINSTALL "\n"
echo -e "\n"
# If some dependencies are missing, asks if user wants to install
if [ "$PKGSTOINSTALL" != "" ]; then
	echo -n "Some dependencies are missing. Want to install them? (y/n): "
	read SURE
	# If user want to install missing dependencies
	if [[ $SURE = "Y" || $SURE = "y" || $SURE = "" ]]; then
		# Debian, Ubuntu and derivatives (with apt-get)
		if hash apt-get &> /dev/null; then
			#apt-get install $PKGSTOINSTALL zlib1g-dev python-numpy libpng-dev java-1.7.0-openjdk-dev libreadline-dev readline-devel bzip2-devel gcc-c++ python-dev python-bzutils libbz2-dev ncurses-dev libncurses5-dev ncurses libncursesw5-dev zlib-dev libbz2-1.0 libbz2-ocaml libbz2-ocaml-dev perl-CPAN
			apt-get install -f $PKGSTOINSTALL java-1.7.0-openjdk-dev libpng-dev build-essential dos2unix gfortran libgd2-xpm-dev libice-dev libncurses5-dev libpdf-api2-perl libpdf-api2-simple-perl libreadline6-dev libxt-dev python2.7-dev python-dev python-matplotlib python-numpy ttf-mscorefonts-installer unzip zlib1g-dev perl-CPAN libbz2-dev ncurses-dev libncurses5-dev ncurses libncursesw5-dev zlib-dev libbz2-1.0 libbz2-ocaml libbz2-ocaml-dev

		# OpenSuse (with zypper)
		elif hash zypper &> /dev/null; then
			zypper in $PKGSTOINSTALL dos2unix bzip2-devel libpng-dev zlib1g-dev java-1.7.0-openjdk-dev libreadline-dev readline-devel python-devel gcc g++ gcc-c++ python-dev python-bzutils libbz2-dev ncurses-devel libncurses5-dev ncurses libncursesw5-dev zlib-dev libbz2-1.0 libbz2-ocaml libbz2-ocaml-dev perl-CPAN
		# Mandriva (with urpmi)
		elif hash urpmi &> /dev/null; then
			urpmi $PKGSTOINSTALL java-1.7.0-openjdk-dev dos2unix bzip2-devel libpng-dev zlib1g-dev libreadline-dev readline-devel python-devel gcc g++ gcc-c++ python-dev python-bzutils libbz2-dev ncurses-devel libncurses5-dev ncurses libncursesw5-dev zlib-dev libbz2-1.0 libbz2-ocaml libbz2-ocaml-dev perl-CPAN
		# Fedora and CentOS (with yum)
		elif hash yum &> /dev/null; then
			yum install $PKGSTOINSTALL dos2unix python-matplotlib numpy libpdf libXt-devel libICE-devel gd-progs gd-devel gcc-gfortran bzip2-devel libpng-devel java-1.7.0-openjdk-devel zlib1g-devel libreadline-devel readline-devel python-devel gcc g++ gcc-c++ python-devel python-bzutils libbz2-devel ncurses-devel libncurses5-devel ncurses zlib-devel libbz2-1.0 libbz2-ocaml libbz2-ocaml-devel perl-CPAN
			yum groupinstall 'Development Tools' 'Development Libraries'
		# ArchLinux (with pacman)
		elif hash pacman &> /dev/null; then
			pacman -Sy $PKGSTOINSTALL bzip2-devel libpng-dev zlib1g-dev libreadline-dev java-1.7.0-openjdk-dev gcc g++ gcc-c++ python-dev python-bzutils libbz2-dev ncurses-devel libncurses5-dev ncurses libncursesw5-dev zlib-dev libbz2-1.0 libbz2-ocaml libbz2-ocaml-dev perl-CPAN
		# Else, if no package manager has been founded
		else
			# Set $NOPKGMANAGER
			NOPKGMANAGER=TRUE
			echo "ERROR: impossible to found a package manager in your system.
Please, install manually $DEPENDENCIES ."
		fi
		# Check if installation is successful
		if [[ $? -eq 0 && ! -z $NOPKGMANAGER ]] ; then
			echo "All dependencies are satisfied."
		# Else, if installation isn't successful
		else
			if hash java && hash g++ && hash wget 2>/dev/null; then
				clear
				echo -e "Error: impossible to install some missing dependencies.
Please install manually:- $PKGSTOINSTALL."
				echo -e "Note: JAVA and g++ do exist, please ignore from above list."
				echo -e "\n"
				else
				echo -e "Error: impossible to install some missing dependencies.
Please install manually:- $PKGSTOINSTALL."
				echo -e "\n"
			fi
		fi
	# Else, if user don't want to install missing dependencies
		else
		echo "WARNING: Some dependencies may be missing. So, please,
install manually:- $PKGSTOINSTALL"
	echo -e "\n"
	fi
fi

################ End of the Installation ###################
echo -e "\n"
echo "Installation done !"
echo "The CBS-miRSeq system dependencies are installed correctly
now you may proceed with CBS-miRSeq installation script."
echo "bash ./Install.CBS-miRSeq.dependencies.v1.0.sh"
echo "Thank you for using the CBS-miRSeq pipeline."
echo -e "\n"
source ~/.bashrc
exec bash --login


