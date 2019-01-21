#!/usr/bin/Rscript



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

# # USAGE: Rscript ./CBS-miRSeq.Required.Packages.R
# # Author: Kesharwani RK; Email: bioinforupesh2009.au@gmail.com
# # version: 1.0
# # Copyright (c) 2019 Kesharwani RK

cat("#Welcome to CBS-miRSeq R-package installation script...\n")
cat("#Downloading packages..., please have a patience.\n")
cat("#Clean installtion may take a while.\n")

#############################
## Name of entire packages ##
#############################
# load multiple packages by once but first to check if installed ??
Allpkgs= c("wordcloud", "EDASeq", "RColorBrewer","BiocInstaller","edgeR","DESeq2","gplots","VennDiagram","plotrix","ReportingTools","hwriter","lattice","S4Vectors","clusterProfiler", "reactome.db", "ReactomePA", "GO.db", "biomaRt","DOSE","networkD3", "igraph", "magrittr","stringi","KEGGprofile","pathview","plyr","gridExtra","grid","grDevices","XML","rJava","crayon", "HTSFilter", "GOSemSim", "mixOmics")

########################
## packages from cran ##
########################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos="http://cran.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

# # Install and load packages from cran
cranPack=c("wordcloud", "RColorBrewer","gplots","VennDiagram","plotrix","hwriter","lattice","networkD3", "igraph", "magrittr","stringi","gridExtra","plyr","grid","grDevices","XML","rJava", "crayon", "mixOmics")
ipak(cranPack)

########################################
## Install packages from bioconductor ##
########################################

source("http://bioconductor.org/biocLite.R") ## online use
biopak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, dependencies = TRUE, suppressUpdates=F,ask=F,suppressAutoUpdate=F)
  sapply(pkg, require, character.only = TRUE)
}

##load multiple packages by once but first to check if installed ??
bioPack=c("EDASeq", "BiocInstaller","edgeR","DESeq2","ReportingTools","S4Vectors","clusterProfiler", "reactome.db", "ReactomePA", "GO.db", "biomaRt","DOSE","KEGGprofile","pathview", "HTSFilter","GOSemSim")
biopak(bioPack)

## Load packages into session 
cat("#Loading packages...\n")

####################
## Load libraries ##
####################

loaded_pack=suppressWarnings(suppressMessages(lapply(Allpkgs, require,  character.only=T)))


############################################################
## Confirmation if requyired all packages has been loaded ##
############################################################

## checking if all packages are loaded successfully
##cat("Checking if packages are loaded successfully...\n")
fine=grep("TRUE", loaded_pack, perl=TRUE, value=F)
if(length(fine) < length(Allpkgs)) # no. of packages to be loaded
{
  stop("Seems required packages are not loaded/installed;\nidentify the package with the broken dependencies(find something:there is no package called..) and to Install it manually, further launch this script again")
} else {
  message("#Packages successfully loaded.\n")
  message("#Now you can processed with CBS-miRSeq analysis.\n")
  message("#Thank you for using the CBS-miRSeq pipeline.\n")
}

############################## Installation done ###########################



## fix problem of rJava, JAVA_HOME
## http://stackoverflow.com/questions/3311940/r-rjava-package-install-failing
## http://stackoverflow.com/questions/13466777/jni-h-no-such-file-or-directory

## fix problem of XML and openssl
# (for centOS)
## sudo yum install curl curl-devel
## sudo yum -y install libxml2 libxml2-devel 
# (for UBUNTU)
## sudo apt-get update
## sudo apt-get curl curl-devel (for UBUNTU)
## sudo apt-get install libxml2-dev (for UBUNTU)
## sudo apt-get install r-cran-xml
## sudo apt-get install libcurl4-openssl-dev
## sudo apt-get install libssl-dev






