#!/usr/bin/Rscript

## Pipeline script for Diff Expr tags

######################## Diff Expr script ############################################
# Description:
# USAGE: CBS-miRSeq.DE.script_v1.R <miR_expression.txt> <output.dir.path> <analysis sps> <startCol.grpsA> <startCol.grpsB>
# Parameters:
# miR_expression = miR Expression matrix obtained byi ht-seq/featuresCounts(*.txt)
# output_dir_path= Path where results should be produced
# sps= Analysis sps
# name.grpA= condition A
# startCol_grpA= Group First start column number in matrix (usually 2)
# name.grpB= condition B
# startCol_grpB= Group Second start column number in matrix
# LowCountsfeaturesFilterByCPM = pre-filter i.e. filtration of unexpressed or low expressed tags;
# CutoffLowCountsfeaturesFilter = filtration CutoffLowCountsfeaturesFilter either used for LowCountsfeaturesFilterByCPM or Raw counts filter if LowCountsfeaturesFilterByCPM is set to "no"
# RemoveInconsistentFeatures = Remove bad features from the matrix before to perform DE
# ThresholdToRemoveInconsistentFeatures = greter than this cutoff tags will be removed(applied in both model)
# IndependentFiltering = if yes then applied to both model
# date: 12 April 201
# version : v1.0
# Authors: Kesharwani RK email: bioinforupesh2009.au@gmail.com
#######################################################################################

# #######################################################################################
# # # Example:
# miR_expression <-  "/Users/Destiny/Desktop/newDataset2/dataset2Counts.txt"
# output.dir.path <- "/Users/Destiny/Desktop/newDataset2/ResultsDataset2/"
# sps <- "hsa" #Analysis sps
# name.grpA="WT" ## condition A
# startCol.grpA <- as.numeric("2")
# name.grpB="ATM" ## condition B
# startCol.grpB <- as.numeric("5")
# LowCountsfeaturesFilterByCPM <- "yes" ## or "no" # [Recommended = "yes"]
# CutoffLowCountsfeaturesFilter <- as.numeric(1) ## Filter unexpressed or low expressed tags;[recommended CutoffLowCountsfeaturesFilter="5"]
# RemoveInconsistentFeatures <- "no" # or "no" [recommended] ### remove Inconsistent from matrix before to perform DE
# ThresholdToRemoveInconsistentFeatures <- as.numeric(1.5) ##[recommended>1.5] ### less than this threshold tags will be kept
# IndependentFiltering <- "no" ## or "yes"
# plotType <- "pdf" ##or "eps" ## [default]
# pval_Cutoff= as.numeric(0.05)
# padj_Cutoff= as.numeric(0.05)
# log2FC_Cutoff= as.numeric(1.5)

# #######################################################################################

# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
## NOTE 1: if LowCountsfeaturesFilterByCPM is "yes" then please also asign the value for CutoffLowCountsfeaturesFilter () [recommended CutoffLowCountsfeaturesFilter="1"];
## otherwise The Filtering will be performded by the raw counts based on assinged CutoffLowCountsfeaturesFilter ().

## Tip1 : If you do not want any low counts features then assign the CutoffLowCountsfeaturesFilter <- "0".
## Therfore It will only filter non-expressed tags.
## Moreover, assigning the value CutoffLowCountsfeaturesFilter <- "" ; will ignore any prefiltering from your matrix.

## Tip2 : IndependentFiltering i.e. Independent filtering is specially fruitful for edgeR as edgeR has never implemented any independent filtering.

## Tip3 : Enabling the option RemoveInconsistentFeatures="yes"; will eliminate the features with inconsistent expression values from your original matrix prior to DE
## We advised this for only small datasets (3-5 sample in each group)

## NOTE 2 : Please keep in mind, once you enable "yes" for IndependentFiltering; independent filtering will be used for both model (edgeR as well as DESeq2).
## resultant, Native independent filtering used in DESEq2 will no longer apply.
## However, Based on practice, we found minor difference in the results while using
## HTSFilter or genefilter implemented into DESeq2.

# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#rm(list=ls()) # cleaning work space
cat("DE Analysis started....\n")
cat ("date & time:\n")
date()

# input arg
arg <- commandArgs(TRUE)

if(length(arg) < 16)
{
 stop("#ERROR! please supply all inputs.
\nUSAGE: Rscript ./CBS-miRSeq.DE.script_v1.R <miR_expression.txt> <output.dir.path>
<3 letter code of analysis sps>
<condition/group.A> <startCol.grpsA>
<condition/group.B> <startCol.grpsB>
<LowCountsfeaturesFilterByCPM(yes/no)> <CutoffLowCountsfeaturesFilter>
<RemoveInconsistentFeatures(yes/no)> <ThresholdToRemoveInconsistentFeatures>
<IndependentFiltering(yes/no)> <plotType(pdf/eps)>
<pval_Cutoff> <padj_Cutoff> <log2FC_Cutoff>
\n")
}
miR_expression <- arg[1]
output.dir.path <- arg[2]
sps <- arg[3]
name.grpA <- arg[4]
startCol.grpA <- as.numeric(arg[5])
name.grpB <- arg[6]
startCol.grpB <- as.numeric(arg[7])
LowCountsfeaturesFilterByCPM <- arg[8]
CutoffLowCountsfeaturesFilter <- as.numeric(arg[9])
RemoveInconsistentFeatures <- arg[10]
ThresholdToRemoveInconsistentFeatures <- as.numeric(arg[11])
IndependentFiltering <- arg[12]
plotType <- arg[13]
pval_Cutoff <- as.numeric(arg[14])
padj_Cutoff <- as.numeric(arg[15])
log2FC_Cutoff <- as.numeric(arg[16])

# #################################
# ### load packages into session ##
# #################################
#
# #load multiple packages by once but first to check if installed ??
pkgs= c("RColorBrewer","BiocInstaller","edgeR","DESeq2","gplots","VennDiagram","plotrix","ReportingTools","hwriter","lattice","HTSFilter","mixOmics","dplyr", "gridBase", "grid", "gridExtra","reshape2","ggplot2", "reshape", "EDASeq", "wordcloud")

# #packages from cran

cranPack=c("RColorBrewer","gplots","VennDiagram","plotrix","hwriter","lattice","mixOmics","dplyr", "gridBase", "grid", "gridExtra","reshape2","ggplot2", "reshape", "wordcloud")

# #Install packages from R (if not already installed)
inst_bio2 <- cranPack %in% installed.packages()
if(length(cranPack[!inst_bio2]) > 0)
{
 suppressWarnings(suppressMessages(lapply(cranPack[!inst_bio2],install.packages, repos="http://cran.r-project.org")))
 #install.packages(pkgs=x,repos="http://cran.r-project.org")
}
## packages from bioconductor
bioPack=c("BiocInstaller","edgeR","DESeq2","ReportingTools","HTSFilter", "EDASeq")
# Install bioconductor packages (if not already installed)
inst_bio <- bioPack %in% installed.packages()

if(length(bioPack[!inst_bio]) > 0)
{
 source("http://bioconductor.org/biocLite.R") ## connect to biocLit
 suppressWarnings(suppressMessages(lapply(bioPack[!inst_bio], biocLite, dependencies = T,suppressUpdates=T,ask=F,suppressAutoUpdate=T)))
}

## Load packages into session
cat("#loading packages...\n")
loaded_pack=suppressWarnings(suppressMessages(lapply(pkgs, require,  character.only=T)))

## checking if all packages are loaded successfully
fine=grep("TRUE", loaded_pack, perl=TRUE, value=F)
if(length(fine) < length(pkgs)) # no. of packages to be loaded
{
 stop("Seems required packages are not loaded into session;\nidentify the package with the broken dependencies or try a different versions to see if that works.")
} else {
 cat("Packages successfully loaded.\n")
}

################################
### 1. load count data matrix ##
########################## ####

cat("#loading digital counts matrix for Differential expression analysis.....\n")
dataset <- read.table(miR_expression, header= TRUE, row.names = 1)
cat ("#No. of Total Tags Before filter:\n")
dataset <- dataset
dim(dataset)

##########################
#### 2. Pre-Filtering ####
##########################

pre_filter <- function (DF, CutoffLowCountsfeaturesFilter = as.numeric(1)) {
  if (LowCountsfeaturesFilterByCPM=="yes") {
    ###Method1  <- filtering by using default cpm counts
    cat ("#Filtering low counts tags by CPM\n")
    cat("#Keeping tags that have CutoffLowCountsfeaturesFilter (Counts per million) in atleast half (50%) of the Total sample length..\n")
    cat("#filtering CutoffLowCountsfeaturesFilter:", CutoffLowCountsfeaturesFilter)
    keep <- rowSums(cpm(DF)> as.numeric(CutoffLowCountsfeaturesFilter)) >= ncol(DF)/2
    data_filt <- DF[keep, ]
  } else {
    ###Method2  <- filtering by using raw counts
    cat ("#Filtering low counts tags by raw counts\n")
    cat("#Keeping tags that have CutoffLowCountsfeaturesFilter in atleast half (50%) of the Total sample length..\n")
    keep <- rowSums(DF> as.numeric(CutoffLowCountsfeaturesFilter)) >=  ncol(DF)/2
    data_filt <- DF[keep, ]
  }
}

## pre-filtering
data_filt <- pre_filter(DF = dataset, CutoffLowCountsfeaturesFilter = CutoffLowCountsfeaturesFilter)
dim(data_filt)


## just for dot plot what was filter from here to cv
data_filt2 <- pre_filter(DF = dataset, CutoffLowCountsfeaturesFilter = CutoffLowCountsfeaturesFilter)
## Dimension of Filtered data
cat ("\n#No. of Tags after filter:\n")
dim(data_filt)
## Total sample
num.samples <- ncol(data_filt)
cat("#Total no. of samples:", num.samples, "\n")
# no. of samples in each condition (group)
num.group.A <- startCol.grpB - startCol.grpA
cat("#No. of samples in group A:", num.group.A, "\n")
num.group.B <- num.samples - startCol.grpB + 2
cat("#No. of samples in group B:", num.group.B, "\n")
##### Inconsistent cfiltering before analysis ####
## function to compute CV
co.va <- function(myDF) {
  sd(myDF)/mean(myDF)
}
RemovebadTags <- function (DF, threshold = as.numeric(1.5), matrix.left=c(keep, out)) {
  group.A <- as.data.frame(DF[,1:num.group.A])
  group.B <- as.data.frame(DF[,(num.group.A+1):num.samples])
  group.A$coVar.A <- apply(cpm(group.A)[,c(colnames(group.A))],1,co.va)
  group.B$coVar.B<-apply(cpm(group.B)[,c(colnames(group.B))],1,co.va)
  all.grp <- merge(as.data.frame(group.A), as.data.frame(group.B), by.x="row.names", by.y="row.names")
  if (matrix.left=="keep") {
    data_filt.keep <- all.grp[(all.grp$coVar.A < ThresholdToRemoveInconsistentFeatures & all.grp$coVar.B < ThresholdToRemoveInconsistentFeatures), ]
  } else
    data_filt.left <- all.grp[(all.grp$coVar.A > ThresholdToRemoveInconsistentFeatures | all.grp$coVar.B > ThresholdToRemoveInconsistentFeatures), ]
}

##################################
#### 3. remove  hypervariants ####
##################################

## Inconsistent filtration
if(RemoveInconsistentFeatures=="yes") {
  data_filt <- RemovebadTags(DF = data_filt, threshold = ThresholdToRemoveInconsistentFeatures, matrix.left = "keep")
  data_filt <- data_filt[!(rowSums(is.na(data_filt))),]
  row.names(data_filt) <- data_filt$Row.names
  data_filt$Row.names <- NULL
  data_filt$coVar.A <- NULL
  data_filt$coVar.B <- NULL
  cat("\n#kept tags that have cv lower than threshold:", ThresholdToRemoveInconsistentFeatures ,"in each class/group..\n")
  dim(data_filt)
} else {
  cat("\n#Computation of coef.variation is not set.")
  cat("\n#Keeping dimension of Matrix is the same as pre-filter.\n")
  dim(data_filt)
}

## extract cv filtered tags
if(RemoveInconsistentFeatures=="yes") {
  cv.out <- RemovebadTags(DF = data_filt2, threshold = ThresholdToRemoveInconsistentFeatures, matrix.left = "out")
  cv.out <- cv.out[!(rowSums(is.na(cv.out))),]
  row.names(cv.out) <- cv.out$Row.names
  cv.out$Row.names <- NULL
  cat("#No. of tags filtered by cv:")
  cat(length(row.names(cv.out)))
  cat("\n")
}

## Extract the tags which have been filtered via CV and expression vs group plot for each removed tags

## Dot plot function for Inconsistent
dot_plot <- function (df, output_name, output.dir.path) {
  plot_list <- list()
  for (i in 1:nrow(df)){
    ### only for one sample
    singleGene <- melt(df[i,], id="Gene")
    colnames(singleGene)[2] <- "Sample"
    colnames(singleGene)[3] <- "Expression"
    singleGene$group <- c(rep(name.grpA, num.group.A), rep(name.grpB, num.group.B))
    lable <- unique(singleGene$Gene)
    plot_list [[i]] <- ggplot(data=singleGene, aes(x=group, y=Expression, shape = Sample, color=Sample)) + scale_shape_manual(values=seq(0,num.samples)) + geom_point(size = 3) + ggtitle(lable)
  }
  mulitiPlot <- marrangeGrob(plot_list, nrow=1, ncol=2)
  ggsave(paste(output.dir.path,"/",output_name,".pdf",sep=""), plot = mulitiPlot, width = 11, height = 8.5, units = "in", dpi = 600)
  dev.off()
}

## extract the tags which have been filtered via CV
if(RemoveInconsistentFeatures=="yes") {
  tagsCV <- rownames(cv.out)
  filtered.tags.by.cv <- data_filt2[rownames(data_filt2) %in% tagsCV, ]
  filtered.tags.by.cv <- as.data.frame(filtered.tags.by.cv)
  if(length(rownames(filtered.tags.by.cv)) > 1) {
    cat("\nwritting an expression table of inconsistent tags....\n")
    write.csv(as.data.frame(cv.out), paste(output.dir.path,"/","filtered.inconsistent.tags.by.cv.csv",sep=""),quote=F,row.names=T)
    cat("Plotting Probable inconsistent tags....\n")
    filtered.tags.by.cv$Gene <- row.names(filtered.tags.by.cv)
    dot_plot(df = filtered.tags.by.cv, output_name = "filtered.inconsistent.tags.by.cv", output.dir.path = output.dir.path)
    cat("done.")
  } else {
    cat("NO Inconsistent tags FOUND at THRESHOLD:", ThresholdToRemoveInconsistentFeatures)
    NoteNoCV <- c("NO Inconsistent tags found at threshold For Remove Inconsistent tags From Matrix:",ThresholdToRemoveInconsistentFeatures,"Seems good expression data.")
    write.table(NoteNoCV, paste(output.dir.path,"/",sps,"_NO.Inconsistent.FOUND.txt",sep=""),quote=F,row.names=F, col.names = F)
  }
}

# # Identifying differentially expressed tags

##################################
### 4. edgeR #TMM normalization ##
##################################

## DGEList object
cat ("\n")
cat("#edgeR Diff expression analysis is running.....\n")

group = as.factor(c(rep(name.grpA, num.group.A), rep(name.grpB, num.group.B)))

###Construction of Object DGEList
DGE.object = DGEList(data_filt, group = group)

### TMM normalization
TMM_norm <- calcNormFactors(DGE.object, method="TMM")
## Note that the “size factor” from DESeq is not equal to the “norm factor” in the edgeR. In edgeR, the library size and additional normalization scaling factors are separated. See the two different columns in the $samples element of the 'DGEList' object above. In all the downstream code, the lib.size and norm.factors are multiplied together to act as the effective library size; this (product) would be similar to DESeq's size factor.
## Estimating normalized absolute expression from their scaling factors
TMM_norm.ScaleFactors <- DGE.object$samples$lib.size * TMM_norm$samples$norm.factors
## To get extact TMM normalized expression
TMM_norm.exp <- (t(t(TMM_norm$counts)/TMM_norm.ScaleFactors) * mean(TMM_norm.ScaleFactors))

# experimental design
design <- model.matrix(~group) ### 'Intercept=control="c"'
cat("#Estimating the Negative Binomial Dispersion..\n")
y <- estimateGLMCommonDisp(TMM_norm, design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

cat("#Running GLM Model for DE analysis...\n")
fit <- glmFit(y,design)
# fit <- glmFit(y,dispersion=y$tagwise.dispersion, design)
lrt <- glmLRT(fit,coef=2)

###############################
## 5. Implementing HTSFilter ##
###############################

if (IndependentFiltering=="yes") {
  cat("#Running HTSFilter for edgeR analysis.....\n")
  lrt <- HTSFilter(lrt, DGEGLM=fit, plot=FALSE)$filteredData
  #dim(lrt)
  # Fisher's exact test for all tags (after HTSFilter)
  Toplrt <- topTags(lrt, n=nrow(lrt$table))
  ## writing HTSFilt table
  edgeR_Results <- merge(as.data.frame(TMM_norm.exp), as.data.frame(Toplrt), by.x="row.names", by.y="row.names",sort=FALSE)
  names(edgeR_Results)[1] <- "Features"
  #order all counts by FDR
  edgeR_Results <- edgeR_Results[order(edgeR_Results$FDR),]
  write.csv(as.data.frame(edgeR_Results), paste(output.dir.path,"/",sps,"_edgeR_HTSFilt.TMM.csv",sep=""),quote=F,row.names=F)
  cat("Done.\n")
} else {
  # Fisher's exact test for all tags
  Toplrt <- topTags(lrt, n = nrow(y), adjust.method = "BH")
  # writing edgeR result in a data frame
  ## Merge with TMM normalized count data
  cat("#Writing edgeR results....\n")
  cat("#Merging TMM normalized counts into Results...")
  edgeR_Results <- merge(as.data.frame(TMM_norm.exp), as.data.frame(Toplrt), by.x="row.names", by.y="row.names",sort=FALSE)
  names(edgeR_Results)[1] <- "Features"
  #order all counts by FDR
  edgeR_Results <- edgeR_Results[order(edgeR_Results$FDR),]
  write.csv(as.data.frame(edgeR_Results), paste(output.dir.path,"/",sps,"_edgeR_results_TMM.csv",sep=""),quote=F,row.names=F)
  cat("\n#Done.\n")
}

################## edgeR LRT HTML Report ####################
cat ("#Creating HTML for edgeR results..\n")
rep.theme <- reporting.theme()
## Change symbol colors in plots
rep.theme$superpose.symbol$col <- c("blue", "red")
rep.theme$superpose.symbol$fill <- c("blue", "red")
lattice.options(default.theme = rep.theme)

deReport <- HTMLReport(shortName = 'CBS-miRSeq.Analysis.with.edgeR',
                       title = 'DE analysis using edgeR(LRT)',
                       basePath=output.dir.path,reportDirectory = "./HTMLRports")

addRowLink <- function(df, ...)
{
 df$IDs <- hwrite(as.character(df$IDs),
                  link = paste0("http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=",as.character(df$IDs)), table = FALSE)
 return(df)
}

publish(lrt, deReport, pvalueCutoff=1, lfc=0, n=nrow(lrt),
        annotation.db=" ", countTable=TMM_norm.exp,
        conditions=c(rep("groupA", num.group.A), rep("groupB", num.group.B)),
        .modifyDF= list(modifyReportDF, addRowLink),
        make.plots = TRUE)
finish(deReport)
cat("\n")

#################################
### 6. Plotting edgeR results ###
#################################

###Total of DE at given CutoffLowCountsfeaturesFilter (FDR) ## up and down regulation check
##cat ("#Tags are differentially expressed at FDR:", padj_Cutoff,"\n")
##t(summary(dt <- decideTestsDGE(lrt, p = as.numeric(padj_Cutoff))))### t=transverse
dt <- decideTestsDGE(lrt, p = as.numeric(padj_Cutoff), adjust.method="BH") # at FDR
# # pick-out which tags are DE
isDE <- as.logical(dt)
DEnames <- rownames (y) [isDE]

########### Extract the normalized expressions by cpm with log=TRUE
## edgeR
edgeRCPM.log <- cpm(data_filt, normalized.lib.sizes=T, log=TRUE)
write.csv(as.data.frame(edgeRCPM.log), paste(output.dir.path,"/",sps,"_edgeRCPM.log2.csv",sep=""),quote=F)
Updown.exp <- edgeRCPM.log[DEnames, ]
#############################################################

### function to plot; defined by user choice

plotingtype <- function (plotType, analysis) {
  if ( plotType=="eps" ) {
    ## plot png
    setEPS()
    postscript(paste(output.dir.path,"/",sps,"_",analysis,"_analysis_plots%03d.eps",sep=""), onefile = F)
    #dev.off()
  } else if ( plotType == "pdf" |  plotType == "" ) {
    # plot pdf
    pdf(paste(output.dir.path,"/",sps,"_",analysis,"_analysis_plots.pdf",sep=""), width = 14, height = 12)
    #dev.off()
  } else {
    stop("Please specify a plot type !!")
  }
}

### function to plot volcano for both test
PlotVolcanoForBoth<- function (DE_DF, pval_Cutoff, log2FC_Cutoff, padj_Cutoff, main="", legpos="top") {
  ## Valcano plots
  with(as.data.frame(DE_DF), plot(logFC, -log10(PValue), col=ifelse((PValue < as.numeric(pval_Cutoff)), "red", "dimgray"), pch=19, main=main))
  ## add darkturquoise color (PValue < 0.05 & abs(logFC)>= 1))
  with(subset(as.data.frame(DE_DF), PValue < as.numeric(pval_Cutoff) & abs(logFC)>= as.numeric(log2FC_Cutoff)), points(logFC, -log10(PValue), pch=19, col="darkturquoise"))
  ## add navy color (FDR < 0.05 & abs(logFC)>= 1)
  with(subset(as.data.frame(DE_DF), FDR < as.numeric(padj_Cutoff) & abs(logFC)>= as.numeric(log2FC_Cutoff)), points(logFC, -log10(PValue), pch=19, col="navy"))
  ## logFC line
  abline(v = c(-as.numeric(log2FC_Cutoff), as.numeric(log2FC_Cutoff)), col = "lightsteelblue4", lty = 3)
  ## label points with significant tags
  Sigtags=rownames(subset(as.data.frame(DE_DF), FDR < as.numeric(padj_Cutoff) & abs(logFC)>= as.numeric(log2FC_Cutoff)))
  ## significant tags
  if (length(Sigtags) >0 & length(Sigtags)<= 25) ## lower than this will be labeled
  {
    with(subset(as.data.frame(DE_DF), FDR < as.numeric(padj_Cutoff) & abs(logFC)>= as.numeric(log2FC_Cutoff)),text(logFC, -log10(PValue), Sigtags, cex=.6, pos=1, col="tan1"))
  }
  ## legends
  mylegend <- c(paste0("PValue > ", pval_Cutoff),
                paste0("PValue < ", pval_Cutoff),
                paste0("PValue < ", pval_Cutoff, " & log2FC", " >= |", log2FC_Cutoff, "|"),
                paste0("Adj.PValue < ", padj_Cutoff, " & log2FC", " >=  |", log2FC_Cutoff, "|"))
  leg.col <- c("dimgray", "red","darkturquoise","navy")
  legend(legpos, legend=rev(mylegend), pch = 19,col = rev(leg.col),cex=.6)
}

cat("#plotting edgeR results...\n")
plotingtype(plotType = plotType, analysis = "edgeR")
## color
col = c(rep("mediumseagreen", times=num.group.A),rep("indianred2", times=num.group.B))
opar <- par()
par(mfrow=c(2,2))
cat("Relative Log Expression (RLE) plot..\n")
suppressWarnings(EDASeq::plotRLE(edgeRCPM.log, col=col, isLog=T, outline=FALSE, ylab="Expression(logCPM)", main="Filtered RawCounts"))
suppressWarnings(EDASeq::plotRLE(TMM_norm.exp, col=col, isLog=F, outline=FALSE, main="TMM Normalization after filter", ylab="TMM_Expression(log2)"))
cat("#Plotting PCA...\n")
suppressWarnings(EDASeq::plotPCA(object = edgeRCPM.log, isLog=T, col=col,cex=1.2))
suppressWarnings(EDASeq::plotPCA(object = TMM_norm.exp, isLog=F, col=col,cex=1.2))
suppressWarnings(par(opar))
# # pvalue hist plot
cat("#Plotting Pvalue Histogram plots...\n")
hist(lrt$table$PValue, breaks=20, col= "grey", border="slateblue", main="", xlab= "edgeR: pvalue distribution")

## Plot volcano for edgeR
PlotVolcanoForBoth(DE_DF = Toplrt$table, pval_Cutoff = pval_Cutoff, log2FC_Cutoff = log2FC_Cutoff, padj_Cutoff = padj_Cutoff)

## plotting heatmap top FDR tags
cat("#Plotting heatmap for Updown tags ...\n")
m=paste0("TopTags at FDR < ", padj_Cutoff)
colors <- c("lemonchiffon3", "gainsboro")
col2=c(rep("lightsteelblue3", times=num.group.A),rep("peachpuff2", times=num.group.B))
if (length(rownames(Updown.exp)) < 50 & length(rownames(Updown.exp)) > 8) {
  ## heatmap
  heatmap.2(Updown.exp, col = greenred(75), scale = "row", ColSideColors = col2, key = T,
            symkey = FALSE, density.info = "none", trace = "none", cexRow = 1,
            cexCol=1,margin = c(8,11), main = m, dendrogram = "both",
            srtCol=45, distfun = function(x) as.dist(1 - cor(t(x), method = "spearman")),
            hclustfun = function(x) hclust(x, method = "average"))
  mtext("[dist: 1-corr(spearman), method: average]", side = 3, adj = 0.6, col = "blue", line=1)
} else if (length(rownames(Updown.exp)) >= 50) {
  ## heatmap
  heatmap.2(Updown.exp[1:50,], col = greenred(75), scale = "row", ColSideColors = col2, key = T,
            symkey = FALSE, density.info = "none", trace = "none", cexRow = 1,
            cexCol=1,margin = c(8,11), main = m, dendrogram = "both",
            srtCol=45, distfun = function(x) as.dist(1 - cor(t(x), method = "spearman")),
            hclustfun = function(x) hclust(x, method = "average"))
  mtext("[dist: 1-corr(spearman), method: average]", side = 3, adj = 0.6, col = "blue", line=1)
} else {
  cat("#NO Heatmap plotted: No. of Top tags","( at FDR<",padj_Cutoff,")","are less than 8\n")
}

graphics.off()

cat("#edgeR analysis done.\n")
cat("\n")

################ END of edgeR analysis ####################

#####################################################
## 7. Differential Expression analysis with DESeq2 ##
#####################################################

cat("#DESeq2 diff expression analysis is running.....\n")
## experiment design
desq.design <- data.frame(row.names=colnames(data_filt), condition = group)

# count data matrix for deseq2
dds <- DESeqDataSetFromMatrix(countData = data_filt, colData=desq.design, design = ~ condition)

cat("#DE analysis based on the Negative Binomial:", "test = LRT\n")
## DESeq function with GLM fitting
ddsLRT = DESeq(dds, test ="LRT", reduced = ~1)
##ddsLRT = DESeq(dds, test ="Wald", betaPrior=T) ## non shrinkage fold changes
#########################################
## 8. HTSfilter / genefilter by DESeq2 ##
#########################################

## performing HTSFilter for DESeq analysis
if (IndependentFiltering=="yes") {
  cat("#Running HTSFilter for DESeq analysis...\n")
  HTSFilter.ddsLRT <- HTSFilter(ddsLRT, plot=FALSE)$filteredData
  cat("#HTSFiltering done.\n")
  cat("#Extract results for HTSFiltered tags..\n")
  # # results out of HTSFilter
  resLRT=results(HTSFilter.ddsLRT, cooksCutoff=T, independentFiltering=FALSE)
  ## writing DESeq2 analysis results
  cat("#Writing deseq2 analysis results....\n")
  ## deseq Normalized values
  NormalizedCountData <- counts(HTSFilter.ddsLRT, normalized=TRUE)
  DESeq2.results <- merge(as.data.frame(NormalizedCountData),as.data.frame(resLRT), by.x='row.names', by.y='row.names', sort=FALSE)
  names(DESeq2.results)[1] <- 'Tags'
  #order all counts by FDR
  DESeq2.results <- DESeq2.results[order(DESeq2.results$padj),]
  write.csv(DESeq2.results, paste(output.dir.path,"/",sps,"_Deseq2_HTSFilter.csv",sep=""),quote=F,row.names=F)
  ## Using variance stabilization to transform data
  vsd <- varianceStabilizingTransformation(HTSFilter.ddsLRT, blind=T)
  # save Transformed matrix
  write.csv(as.data.frame(assay(vsd)), paste(output.dir.path,"/",sps,"_DESeq2_vst_transformed_counts.csv",sep=""),quote=F)
} else {
  # # results out without HTSFilter
  ## If HTSFilter off then enable native "independentFiltering" from DESeq2
  cat("#Extract results from native deseq analysis..\n")
  resLRT=results(ddsLRT, cooksCutoff=T, independentFiltering=TRUE)
  ## deseq Normalized values
  cat("#Writing deseq2 analysis results....\n")
  NormalizedCountData <- counts(ddsLRT, normalized=TRUE)
  DESeq2.results <- merge(as.data.frame(NormalizedCountData), as.data.frame(resLRT),by.x='row.names', by.y='row.names', sort=FALSE)
  names(DESeq2.results)[1] <- 'Tags'
  #order all counts by FDR
  DESeq2.results <- DESeq2.results[order(DESeq2.results$padj),]
  write.csv(DESeq2.results, paste(output.dir.path,"/",sps,"_DESeq2_results.csv",sep=""),quote=F,row.names=F)
  ## Using variance stabilization to transform data
  vsd <- varianceStabilizingTransformation(ddsLRT, blind=T)
  # save Transformed matrix
  write.csv(as.data.frame(assay(vsd)), paste(output.dir.path,"/",sps,"_DESeq2_vst_transformed_counts.csv",sep=""),quote=F)
  cat("#Done.\n")
}

############### DESEq2 HTML report #####################
cat ("#Creating HTML for DESeq2 results..\n")
des2Report <- HTMLReport(shortName = 'CBS-miRseq.Analysis.with.DESeq2',
                         title = 'DE analysis using DESeq2(LRT)',
                         basePath=output.dir.path,reportDirectory = "./HTMLRports")

addRowLink <- function(df, ...)
{
 df$ID <- hwrite(as.character(df$ID),
                 link = paste0("http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=",as.character(df$ID)), table = FALSE)
 return(df)
}

publish(resLRT, des2Report, pvalueCutoff=1, lfc=0, n=nrow(resLRT),
        annotation.db=" ", DataSet = NormalizedCountData,
        factor = dds$condition,
        .modifyDF= list(modifyReportDF, addRowLink),
        make.plots = TRUE)

finish(des2Report)

################################
## 9. Plotting DESeq2 results ##
################################

# Heatmap of top-hit tags in deseq2 (padj given by user)
## sorted by padj
resLRT <- resLRT[order(resLRT$padj), ]
## extraction of DE at value < padj_Cutoff
toptags=row.names(resLRT)[resLRT$padj < padj_Cutoff]
#is.na(miRNAs.DESeq2)
toptags=toptags[ !is.na(toptags) ]
#sum(na.omit(resLRT$padj < 0.05)) ### pvalue
## Transform raw counts into normalized values (Transformed matrix)
toptags.vsd <- as.data.frame(assay(vsd))[toptags, ]
cat ("#plotting analysis of DESeq2 results...\n")
plotingtype(plotType = plotType, analysis = "DESeq2")
opar <- par()
par(mfrow=c(2,2))
## Raw counts box plots
cat("Relative Log Expression (RLE) plot..\n")
suppressWarnings(EDASeq::plotRLE(edgeRCPM.log, col=col,isLog=T,outline=FALSE,ylab="Expression(logCPM)", main="Filtered RawCounts"))
suppressWarnings(EDASeq::plotRLE(NormalizedCountData, col=col,isLog=F,outline=FALSE, main="DESeq2 Normalization after filter", ylab="DESeq2 Normalized Expression(log2)"))
cat("#Plotting PCA...\n")
suppressWarnings(EDASeq::plotPCA(object = edgeRCPM.log, isLog=T, col=col,cex=1.2))
suppressWarnings(EDASeq::plotPCA(object = NormalizedCountData, isLog=F, col=col,cex=1.2))
suppressWarnings(par(opar))
##pvalue hist plot
## DE histogram at raw pvlaue < 0.05
hist(resLRT$pvalue, breaks=20, col="darkgrey", border="slateblue", main="", xlab= "DESeq2:P-value distribution")

## volcano plot for deseq2 results
DESeqOut <- as.data.frame(resLRT)
## rename column name
colnames(DESeqOut) <- c("baseMean", "logFC", "lfcSE", "stat", "PValue", "FDR")
## plot
PlotVolcanoForBoth(DE_DF = DESeqOut, pval_Cutoff = pval_Cutoff, log2FC_Cutoff = log2FC_Cutoff, padj_Cutoff = padj_Cutoff)

## Clustering and Heatmap
colors <- c("lemonchiffon3", "gainsboro")
col2=c(rep("lightsteelblue3", times=num.group.A),rep("peachpuff2", times=num.group.B))
m=paste0("topTags at FDR < ", padj_Cutoff)
if (length(rownames(toptags.vsd)) < 50 & length(rownames(toptags.vsd)) > 8) {
  ## heatmap
  heatmap.2(as.matrix(toptags.vsd), col = greenred(75), scale = "row", ColSideColors = col2, key = T,
            symkey = FALSE, density.info = "none", trace = "none", cexRow = 1,
            cexCol=1,margin = c(8,11), main = m, dendrogram = "both",
            srtCol=45, distfun = function(x) as.dist(1 - cor(t(x), method = "spearman")),
            hclustfun = function(x) hclust(x, method = "average"))
  mtext("[dist: 1-corr(spearman), method: average]", side = 3, adj = 0.6, col = "blue", line = 1)
} else if (length(rownames(Updown.exp)) >= 50) {
  ## heatmap
  heatmap.2(as.matrix(toptags.vsd[1:50,]), col = greenred(75), scale = "row", ColSideColors = col2, key = T,
            symkey = FALSE, density.info = "none", trace = "none", cexRow = 1,
            cexCol=1,margin = c(8,11), main = m, dendrogram = "both",
            srtCol=45, distfun = function(x) as.dist(1 - cor(t(x), method = "spearman")),
            hclustfun = function(x) hclust(x, method = "average"))
  mtext("[dist: 1-corr(spearman), method: average]", side = 3, adj = 0.6, col = "blue", line = 1)
} else {
  cat("#NO Heatmap plotted: No. of Top tags","( at FDR<",padj_Cutoff,")","are less than 8\n")
}

graphics.off ()

cat("#DESeq2 analysis done.\n")
cat("\n")

################ END of DESeq2 analysis ####################

####################
## 10. Comparison ##
###################

cat ("#Running comparison analysis...\n")

pdf(paste(output.dir.path,"/",sps,"_comparison.plots.pdf",sep=""),width=12,height=12)

#========================================================
##########################################
## 1) intersect by Pavlaue < 0.05 ONLY
##########################################

# compare DE tags by p-value, detected by both models
cat ("#Comparing edgeR and DESeq2 findings at the p-value < 0.05, respectively...\n")
DE.edgeR.PValue <- Toplrt$table$PValue < 0.05
## Number of tags found in edgeR
sum(na.omit(DE.edgeR.PValue))
DE.DESeq2.PValue <- resLRT$pvalue < 0.05
## Number of tags found in edgeR
sum(na.omit(DE.DESeq2.PValue))

## which DE tags are invloved in which model
tags.edgeR.PValue <- row.names(Toplrt$table)[DE.edgeR.PValue]
tags.DESeq2.PValue <- row.names(resLRT)[DE.DESeq2.PValue]
#is.na(tags.DESeq2)
tags.DESeq2.PValue=tags.DESeq2.PValue[ !is.na(tags.DESeq2.PValue) ]

venn(list(edgeR=tags.edgeR.PValue, DESeq2=tags.DESeq2.PValue))
title ("Shared DE tags by p.value (< 0.05)")

#========================================================
#################################################
## 2) intersect by Pvalue < 0.05 and  logF >= |1|
#################################################

# compare DE tags by p-value, detected by both models
cat ("#Comparing edgeR and DESeq2 findings at the p-value < 0.05 and logFC >= |1|, respectively...\n")
DE.edgeR.PValue.LFC1 <- Toplrt$table$PValue < 0.05 & abs(Toplrt$table$logFC) >= 1
## Number of tags found in edgeR
sum(na.omit(DE.edgeR.PValue.LFC1))
DE.DESeq2.PValue.LFC1 <-resLRT$pvalue < 0.05 & abs(resLRT$log2FoldChange) >= 1
## Number of tags found in edgeR
sum(na.omit(DE.DESeq2.PValue.LFC1))

## which DE tags are invloved in which model
tags.edgeR.PValue.LFC1 <- row.names(Toplrt$table)[DE.edgeR.PValue.LFC1]
tags.DESeq2.PValue.LFC1 <- row.names(resLRT)[DE.DESeq2.PValue.LFC1]
#is.na(tags.DESeq2)
tags.DESeq2.PValue.LFC1=tags.DESeq2.PValue.LFC1[ !is.na(tags.DESeq2.PValue.LFC1) ]


venn(list(edgeR=tags.edgeR.PValue.LFC1, DESeq2=tags.DESeq2.PValue.LFC1))
title ("Shared DE tags by p.value(<0.05) and logFC>=|1|")


#========================================================


###################################################
## 3) intersect by Pvalue < 0.05 and  logFC >=|1.5|
###################################################

# compare DE tags by p-value, detected by both models
cat ("#Comparing edgeR and DESeq2 findings at the p-value < 0.05 and logFC >=|1.5|, respectively...\n")
DE.edgeR.PValue.LFC1.5 <- Toplrt$table$PValue < 0.05 & abs(Toplrt$table$logFC) >= 1.5
## Number of tags found in edgeR
sum(na.omit(DE.edgeR.PValue.LFC1.5))
DE.DESeq2.PValue.LFC1.5 <-resLRT$pvalue < 0.05 & abs(resLRT$log2FoldChange) >= 1.5
## Number of tags found in edgeR
sum(na.omit(DE.DESeq2.PValue.LFC1.5))

## which DE tags are invloved in which model
tags.edgeR.PValue.LFC1.5 <- row.names(Toplrt$table)[DE.edgeR.PValue.LFC1.5]
tags.DESeq2.PValue.LFC1.5 <- row.names(resLRT)[DE.DESeq2.PValue.LFC1.5]
#is.na(tags.DESeq2)
tags.DESeq2.PValue.LFC1.5=tags.DESeq2.PValue.LFC1.5[ !is.na(tags.DESeq2.PValue.LFC1.5) ]

venn(list(edgeR=tags.edgeR.PValue.LFC1.5, DESeq2=tags.DESeq2.PValue.LFC1.5))
title ("Shared DE tags by p.value(<0.05) and logFC >= |1.5|")

#========================================================
##########################################
## 4) intersect by FDR < 0.05
##########################################
cat ("#Comparing edgeR and DESeq2 findings at the FDR < 0.05, respectively...\n")
DE.edgeR.FDR <- Toplrt$table$FDR < 0.05
## Number of tags found in edgeR
sum(na.omit(DE.edgeR.FDR))
DE.DESeq2.FDR <-resLRT$padj < 0.05
## Number of tags found in edgeR
sum(na.omit(DE.DESeq2.FDR))

## which DE tags are invloved in which model
tags.edgeR.FDR <- row.names(Toplrt$table)[DE.edgeR.FDR]
tags.DESeq2.FDR <- row.names(resLRT)[DE.DESeq2.FDR]
#is.na(tags.DESeq2)
tags.DESeq2.FDR=tags.DESeq2.FDR[ !is.na(tags.DESeq2.FDR) ]

venn(list(edgeR=tags.edgeR.FDR, DESeq2=tags.DESeq2.FDR))
title ("Shared DE tags by FDR(<0.05)")


#========================================================
############################################
## 5) intersect by FDR < 0.05 and logFC >=|1|
############################################
cat ("#Comparing edgeR and DESeq2 findings at the FDR < 0.05 and logFC >=|1|, respectively...\n")
DE.edgeR.FDR.FC1 <- Toplrt$table$FDR < 0.05 & abs(Toplrt$table$logFC) >= 1
## Number of tags found in edgeR
sum(na.omit(DE.edgeR.FDR.FC1))
DE.DESeq2.FDR.FC1 <- resLRT$padj < 0.05 & abs(resLRT$log2FoldChange) >= 1
## Number of tags found in edgeR
sum(na.omit(DE.DESeq2.FDR.FC1))

## which DE tags are invloved in which model
tags.edgeR.FDR.FC1 <- row.names(Toplrt$table)[DE.edgeR.FDR.FC1]
tags.DESeq2.FDR.FC1 <- row.names(resLRT)[DE.DESeq2.FDR.FC1]
#is.na(tags.DESeq2)
tags.DESeq2.FDR.FC1=tags.DESeq2.FDR.FC1[ !is.na(tags.DESeq2.FDR.FC1) ]

venn(list(edgeR=tags.edgeR.FDR.FC1, DESeq2=tags.DESeq2.FDR.FC1))
title ("#Shared DE tags by FDR(<0.05) and logFC >=|1|")

#========================================================
################################################
## 6) intersect by FDR < 0.05 and  logFC >=|1.5|
################################################
cat ("#Comparing edgeR and DESeq2 findings at the FDR < 0.05 and logFC >=|1.5|, respectively...\n")
DE.edgeR.FDR.FC1.5 <- Toplrt$table$FDR < 0.05 & abs(Toplrt$table$logFC) >= 1.5
## Number of tags found in edgeR
sum(na.omit(DE.edgeR.FDR.FC1.5))

DE.DESeq2.FDR.FC1.5 <- resLRT$padj < 0.05 & abs(resLRT$log2FoldChange) >= 1.5
## Number of tags found in edgeR
sum(na.omit(DE.DESeq2.FDR.FC1.5))

## which DE tags are invloved in which model
tags.edgeR.FDR.FC1.5 <- row.names(Toplrt$table)[DE.edgeR.FDR.FC1.5]
tags.DESeq2.FDR.FC1.5 <- row.names(resLRT)[DE.DESeq2.FDR.FC1.5]
#is.na(tags.DESeq2)
tags.DESeq2.FDR.FC1.5=tags.DESeq2.FDR.FC1.5[ !is.na(tags.DESeq2.FDR.FC1.5) ]

venn(list(edgeR=tags.edgeR.FDR.FC1.5, DESeq2=tags.DESeq2.FDR.FC1.5))
title ("Shared DE tags by FDR(<0.05) and logFC >=|1.5|")

####################################
#### 11. p-value merged results ####
####################################

################ From this point analysis will be performed based on p-value ONLY #########
### The original idea were found and modified accordingly
### http://www.gettinggeneticsdone.com/2012/09/deseq-vs-edger-comparison.html

# tags were only detected by edgeR
de.edgeR.only <- setdiff(tags.edgeR.PValue, tags.DESeq2.PValue)
# tags were only detected by deseq2
de.deseq2.only <- setdiff(tags.DESeq2.PValue,tags.edgeR.PValue)
## common tags in both model
shared.de.tags <- intersect(tags.edgeR.PValue, tags.DESeq2.PValue)

### Extracting logFC, P-values and FDR from both model
## from edgeR
edgeR.out <- as.data.frame(Toplrt$table) ## full results
intersect.stat.edgeR = subset(edgeR.out,row.names(edgeR.out) %in% shared.de.tags)
intersect.stat.edgeR <- intersect.stat.edgeR[,c(1,4,5)]

## from Deseq2
DESeq2.out <- as.data.frame(resLRT) ## full results
intersect.stat.deseq2 = subset(DESeq2.out,row.names(DESeq2.out) %in% shared.de.tags)
intersect.stat.deseq2 <- intersect.stat.deseq2[,c(2,5,6)]

## megring both table into a dataframe
merged.stat <- merge(intersect.stat.edgeR, intersect.stat.deseq2, by.x='row.names', by.y = "row.names", sort= F)

## renaming in order avoid confusions
colnames(merged.stat) <- c("Features","logFC(edgeR)","PValue(edgeR)","FDR(edgeR)","logFC(DESeq2)","Pvalue(DESeq2)","FDR(DESeq2)")

## saving intersect results among all stats
cat ("#Writing an intersect (DE tags found by both methods) results among all stats .....\n")
write.csv(as.data.frame(merged.stat), paste(output.dir.path,"/",sps,"_miR.Intersect.merged.statPvalue0.05.csv",sep=""),quote=F,row.names = F)

##############################
## 12. Full results merged ##
#############################

merged.edgeR.DESeq2 <- merge(edgeR.out, DESeq2.out, by.x ='row.names', by.y = "row.names", sort = FALSE)
colnames(merged.edgeR.DESeq2)[1] <- "Features"

write.csv(as.data.frame(merged.edgeR.DESeq2), paste(output.dir.path,"/",sps,"_merged.edgeR.DESeq2.csv",sep=""),quote=F,row.names = F)

#######################################
## 13. Comparison by the fold changes ##
#######################################

with(merged.edgeR.DESeq2, plot(logFC, log2FoldChange, xlab="edgeR_logFC", ylab="DESeq2_logFC", pch=20, col= "azure4", main="Log2 Fold changes:DESeq2 vs edgeR"))

### edgeR DE tags
with(subset(merged.edgeR.DESeq2, FDR<0.05 & (padj>=0.05 | is.na(padj))), points(logFC,log2FoldChange, pch=20, col="red"))

##deseq DE tags
with(subset(merged.edgeR.DESeq2, padj<0.05 & FDR>=0.05), points(logFC, log2FoldChange, pch=20, col="green"))

## FDR in BOTH (Intersect)
with(subset(merged.edgeR.DESeq2, FDR<0.05 & padj<0.05), points(logFC, log2FoldChange, pch=20, col="blue"))

## legend
legend("topleft", xjust=1, yjust=1, legend=c("FDR<0.05 edgeR", "FDR<0.05 DESeq2", "FDR<0.05(in Both)", "FDR>0.05(in Both)"), pch=19, col=c("red", "green", "blue", "azure4"), bty="n", cex=0.8)
## line
abline(h = c(-1, 1), col = "violet", lty=2)
abline(v = c(-1, 1), col = "orange", lty=2)
graphics.off()

######################
## 14. Word Clouds ##
#####################

cat("Plotting wordclouds of Relative Expression Counts..")
exp.tags.Mean <- rowMeans(data_filt, na.rm = FALSE, dim=1)
exp.tags.Mean <- as.data.frame(exp.tags.Mean)
#miRs <- substring(row.names(exp.tags.Mean), 5)
pdf(paste(output.dir.path,"/",sps,"RelativeExpressionCounts.pdf",sep=""),width=15,height=15, pointsize = 6)
col_clouds <- brewer.pal(8, "Dark2")
wordcloud(row.names(exp.tags.Mean), exp.tags.Mean$exp.tags.Mean, scale = c(8,.8), max.words = Inf, random.order = FALSE, rot.per = .1, colors = col_clouds)
graphics.off()

# # save work space
cat ("#saving Rdata....\n")
save.image(paste(output.dir.path,"/",sps,"_CBS-miRSeq.pipeline.DESeq2vsEdgeR.RData", sep=""))
cat("#Done.\n")
cat("\n#Differiential Analysis has done, please check your Results.\n")
cat("\n#Thank you for using CBS-miRSeq.pipeline.\n")

