#!/usr/bin/Rscript

## Pipeline script for creating Advanced HTML file for expressed known miRNAs with miRBase hyper-link

######################## Diff Expr script ############################################
# Description:
# USAGE: CBS-miRSeq.HTML_v1.R <miRNA_expression.txt> <output.dir.path> <analysis sps>
# Parameters:
# miRNA_expression.txt= miRNA Expression matrix obtained by ht-seq/featuresCounts(*.txt)
# output_dir_path= Path where results should be produced
# sps= Analysis sps
# date: 14 October 2015
# version : v1.0
# Authors: Kesharwani RK email: bioinforupesh2009.au@gmail.com
# Copyright (c) 2019 Kesharwani RK
#######################################################################################

# # ###################################################################################
# example:
# raw.exprs <- "miR.Expression.counts.txt"
# output.dir.path <- "/..../results"
# sps <- "dre"
# # ###################################################################################

#rm(list=ls()) # clean work space if you have anything
cat ("date & time:\n")
date()
# input arg
agr <- commandArgs(TRUE)

if(length(agr) < 3)
{
 stop("please supply all inputs. \nUSAGE: CBS-miRSeq.HTML_v1.R <miRNA_expression.txt> <output.dir.path> <species_3_letter_code(hsa/dre/mmu/rno)>")
}
raw.exprs <- agr[1]
output.dir.path <- agr[2]
sps <- agr[3]

# load multiple packages by once but first to check if installed ??
pkgs= c("ReportingTools","hwriter","lattice","XML")

# Install bioconductor packages (if not already installed)
inst_bio <- pkgs %in% installed.packages()

if(length(pkgs[!inst_bio]) > 0) 
{
 source("http://bioconductor.org/biocLite.R") ## online use
 suppressWarnings(suppressMessages(lapply(pkgs[!inst_bio], biocLite, dependencies = T,suppressUpdates=T,ask=F,suppressAutoUpdate=T)))
}
## Load packages into session 
cat("loading packages...\n")
loaded_pack=suppressWarnings(suppressMessages(lapply(pkgs, require,  character.only=T)))

## checking if all packages are loaded successfully
fine=grep("TRUE", loaded_pack, perl=TRUE, value=F)
if(length(fine) < length(pkgs)) # no. of packages to be loaded
{
 stop("Seems required packages are not loaded into session;\nidentify the package with the broken dependencies or try a different versions to see if that works.")
} else {
 cat("Packages successfully loaded.\n")
}

####################### Analysis started ###########################

miRExp <- read.table(raw.exprs,header= T)
#dim(miRExp)
miRExp_Report <- HTMLReport(shortName = 'CBS-miRseq_miR_Expr',
                            title = 'Digital known miR Expression counts',
                            basePath=output.dir.path,reportDirectory = "./miRExpressionHTMLRports")

addRowLink <- function(df, ...)
{
 df$miRNAs <- hwrite(as.character(df$miRNAs), 
                  link = paste0("http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=",as.character(df$miRNAs)), table = FALSE)
 return(df)
}

publish(miRExp, miRExp_Report, .modifyDF= list(modifyReportDF, addRowLink),make.plots = FALSE)

finish(miRExp_Report)
cat("\nDifferential expr HTML creation has successfully done, go to your output.dir.path to check results.\n\n")






