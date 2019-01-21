#!/usr/bin/Rscript

## Pipeline script for gene ontology (GO), Path and  Network analysis

######################## gene ontology analysis #######################################
# # Description:
# # option should be equal to wrapper shell script
# # Parameters:
# # intersect.genes.uniq.pair = "intersect.genes.uniq.pair.txt" obtained from TargetProdcition
# # by CBS-miRSeq module
# # output_dir_path= Path where results should be produced
# # sps = your species name (3 letter code ONLY)
# # organism <- ## full name of the analysis sps(organism)
# # OrgDb <- ## annotation PACKAGE for your analysis organism
# # entrezID = your gene if entrezID (recommended) or gene symbol
# # GeneID = GeneID should be one of: SYMBOL,REFSEQ,ENSEMBL (## case sensitive)
# # targetHub = "small" (recommended) or "big" or "FullNetwork" ## case sensitive
# # pathways <- "reactome" or "kegg" ## case sensitive
# # download_latest_kegg <- "yes" or "no" if internet working on your working terminal ## case sensitive
# # plotType <- "pdf" or "ps" type in order to plot GO
# # pvalueCutoff <- #the minimum p value for enriched GO and pathways
# # qvalueCutoff <- #the minimum adjusted p value for enriched GO and pathways
# # Updated date : 11 dec 2016
# # Copyright (c) 2019 Kesharwani RK
# #######################################################################################
##cat("\014")
rm(list=ls()) # clean work space if you have anything
cat("#GO enrichment, Pathways and Network analysis started....\n")
cat ("date..\n")
date()
# # input arg
agr <- commandArgs(TRUE)

if(length(agr) < 14)
{
  stop("Incorrect number of arguments ! \nUSAGE:
Rscript ./CBS-miRSeq.GO_analysis.R <miRNA-Target-pair.txt> <output_dir> <organism(NAME OF THE species)> <species 3 letter code> <OrgDb(annotationDB like org.Hs.eg.db)> <entrezID(YES/NO)> <targetHub(small/big/FullNetwork)> <pathways(reactome/kegg)> <download_latest_keg(yes/no)> <plotType(pdf/ps)> <Ontology(ANY/DO)> <GeneID(SYMBOL/REFSEQ/ENSEMBL)> <pvalueCutoff> <qvalueCutoff>\n")
}

## assign inputs
intersect.genes.uniq.pair <- agr[1]
output.dir.path <- agr[2]
organism <- agr[3]
sps <- agr[4]
OrgDb <- agr[5]
entrezID <- agr[6]
targetHub <- agr[7]
pathways <- agr[8]
download_latest_kegg <- agr[9]
plotType <- agr[10]
Ontology <- agr[11]
GeneID <- agr[12]
pvalueCutoff <- agr[13]
qvalueCutoff <- agr[14]
#
# # # # # # # # # # ########################################################################
# # # example:
# intersect.genes.uniq.pair <- "/Users/Destiny/Desktop/WorkingGO/PredictedPairs.txt"
# output.dir.path <- "/Users/Destiny/Desktop/WorkingGO/outGO/"
# sps <- "mmu" ## three letter code of organism
# organism <- "mouse" ## full name of your analysis sps(organism)## case sensitive[small letter only]
# OrgDb <- "org.Mm.eg.db" ## annotation PACKAGE for your analysis organism
# entrezID <- "no" ## YES or NO # case sensitive
# targetHub <- "small" ## should be "small" or "big" or "FullNetwork" ## case sensitive
# pathways <- "kegg" ## or "reactome" or "kegg"
# download_latest_kegg <- "yes" ## should be "yes" or "no" (Needs internet connectivity)
# plotType <- "pdf" ## pdf or ps
# Ontology <- "any" ## "DO" or "ANY" ; ANY (i.e. BP,MF,CC) ==> means except DO (disease ontology)
# GeneID <- "SYMBOL" ##ID should be one of: SYMBOL,REFSEQ,ENSEMBL (## case sensitive)
# pvalueCutoff <- 0.05 #the minimum p value for enriched GO and pathways
# qvalueCutoff <- 0.1  #the minimum adjusted p value for enriched Gene
# ######################################################################
##Install annotation packages (if not already installed)
anno <-c("GO.db", OrgDb)
annopack <- anno %in% installed.packages()
if(length(anno[!annopack]) > 0)
{
  source("http://bioconductor.org/biocLite.R") ## on-line use
  suppressWarnings(suppressMessages(lapply(anno[!annopack], biocLite, dependencies = T,suppressUpdates=F,ask=F,suppressAutoUpdate=F)))
}

annoLoad <- suppressWarnings(suppressMessages(lapply(anno , require,  character.only=T)))
## checking if annotation packages are loaded successfully
fine2=grep("TRUE", annoLoad, perl=TRUE, value=F)
if(length(fine2) < length(anno)) # no. of annotation packages to be loaded
{
  cat("#Wrong input for : OrgDb or organism or sps !!\n")
  cat("#Please recheck and run the script again.\n")
  stop("#OR, Given organism is unable to download, please install manually from http://www.bioconductor.org/packages/3.0/data/annotation/ .\n")
}

# # load multiple packages by once but first to check if installed ??
pkgs= c("S4Vectors","clusterProfiler","ReactomePA", "GO.db", "biomaRt","DOSE","KEGGprofile","pathview","DO.db")
pkgs2=c("networkD3", "igraph", "magrittr","stringi","gridExtra","plyr","grid","grDevices")
# #Install bioconductor packages (if not already installed)
inst_bio <- pkgs %in% installed.packages()
if(length(pkgs[!inst_bio]) > 0)
{
  source("http://bioconductor.org/biocLite.R") ## on-line use
  suppressWarnings(suppressMessages(lapply(pkgs[!inst_bio], biocLite, dependencies =  T,suppressUpdates=F,ask=F,suppressAutoUpdate=F)))
}
# #Install packages from R (if not already installed)
inst_bio2 <- pkgs2 %in% installed.packages()
if(length(pkgs2[!inst_bio2]) > 0)
{
 suppressWarnings(suppressMessages(lapply(pkgs2[!inst_bio2],install.packages, repos="http://cran.r-project.org")))
  install.packages(pkgs=x,repos="http://cran.r-project.org")
}

## all packages
FullPack=c("S4Vectors","clusterProfiler","ReactomePA", "GO.db", "biomaRt","DOSE","networkD3", "igraph", "magrittr","stringi","KEGGprofile","pathview","plyr","gridExtra","grid","grDevices")
## Load packages into session
cat("loading packages...\n")
sapply(FullPack, require,  character.only=T)

loaded_pack=suppressWarnings(suppressMessages(lapply(FullPack, require,  character.only=T)))
## checking if all packages are loaded successfully
fine=grep("TRUE", loaded_pack, perl=TRUE, value=F)
if(length(fine) < length(FullPack)) # no. of packages to be loaded
{
  stop("Seems required packages are not loaded into session;\nidentify the package with the broken dependencies or try a different versions to see if that works.")
} else {
  cat("Packages successfully loaded.\n")
}

########################
## Network Analysis
########################

cat("#GO enrichment, Pathways and Network analysis started....\n")
gene.target = read.table(intersect.genes.uniq.pair, header=T)

## function to counts small and big targets for perticular miRs
huBs <- function(gene.target) {
  res <- table(gene.target$miRNA)
  res <- sort(res, decreasing = TRUE)
  mirs <- names(res)
  resnum <- as.numeric(res)
  freq <- vector()
  for (i in 1:length(mirs)) {
    freq[i] <- paste(gene.target$mRNA[gene.target$miRNA == mirs[i]], collapse = ",")
  }
  res <- data.frame(freq = resnum, names = freq)
  rownames(res) <- mirs
  res <- subset(res, freq >=4)
  return(res)
}

### choose hub
targHub=huBs(gene.target = gene.target)[2]
chooseHub <- function (netPlots) {
  if (targetHub=="small" || targetHub=="Small" || targetHub=="SMALL") {
    smallhub=tail(netPlots,n=1) ## tail + n=1 ; in order to select lower pair
  } else {
    bighub=head(netPlots,n=1) ## tail + n=1 ; in order to select bigger pair
  }
}

# #edges
e=as.vector(targHub$names)
e=chooseHub(e)
## split sep ","
e <-as.vector(strsplit(e, ","))
## unlist object into vector
e <- as.character(unlist(e))
# verteces
v=as.vector(rownames(targHub))
v <- chooseHub(v)
v=rep(v,(length(e)))
## make data frame of target netwrok
userHub=as.data.frame(cbind(v,e))
colnames(userHub) <- colnames(gene.target)

## 3 D html network function
htmlNet <- function (pairs, fileName="") {
  if (targetHub=="FullNetwork" || targetHub=="FULLNETWORK" || targetHub=="fullNetwork") {
    simpleNetwork(pairs, charge = -1000, fontSize = 25,linkColour = "red", nodeColour = "green", nodeClickColour = "darkviolet", textColour = "cornflowerblue", opacity = 1, zoom = TRUE) %>% saveNetwork(file = paste(output.dir.path,"/",sps,"_",fileName,"_3d.network.html",sep=""),selfcontained = F)
  } else {
    simpleNetwork(userHub, charge = -1000, fontSize = 25,linkColour = "red", nodeColour = "green", nodeClickColour = "darkviolet", textColour = "cornflowerblue", opacity = 1, zoom = TRUE) %>% saveNetwork(file = paste(output.dir.path,"/",sps,"_",fileName,"_3d.network.html",sep=""),selfcontained = F)
  }
}

# 3d Plot by networkD3 and generate a html plot
htmlNet(pairs = gene.target, fileName = "miR.mRNA")

## function to create GO enriched output (txt)
print.go.out <- function(goResults, outName) {
  go_rich <- as.data.frame(goResults)
  go_empty <- is.data.frame(go_rich) && nrow(go_rich)==0
  if(go_empty==TRUE) {
    message("Seems some of Your ontology enrichment (BP/MF/CC/pathways) results are empty, please try with lower pvalue and qvalue Cutoff.\n")
  } else {
    ## enriched GO+Description+GenERatio+pvalue
    go <- suppressWarnings(as.data.frame(goResults)[,c(1:3,5:7)])
    #   #order by padj
    go <- go[with(go, order(p.adjust)), ]
    colnames(go)[1] <- "OntologyID"
    #   #rownames(go) <- NULL
    write.table(as.data.frame(go), paste(output.dir.path,"/",sps,outName,"_go.txt",sep=""),quote=F, sep = "\t",row.names = F)
    ## enriched GO and their genes
    ## unlist object
    GO.list <- ldply(goResults@geneSets, data.frame)
    colnames(GO.list) <- c("EnrichedGO","AffliatedGenes")
    write.table(GO.list, paste(output.dir.path,"/",sps,outName,"_GOenriched.txt",sep=""),quote=F,row.names=F,sep = "\t")
  }
}

### function to plot; defined by user choice
### user must provide plot type, there is no defualt option
plotingtype <- function (plotType, main) {
  if ( plotType=="eps" ) {
    ## plot postscript
    setEPS()
    cairo_ps(paste(output.dir.path,"/",sps,"_GOenrich_Plots","_",main,"%03d.eps",sep=""),onefile = F)
    ###postscript(paste(output.dir.path,"/",sps,"_GOenrich_Plots","_",main,"%03d.eps",sep=""), onefile = F)
	#dev.off()
  } else if ( plotType =="pdf" ) {
    # plot pdf
    cairo_pdf(paste(output.dir.path,"/",sps,"_GOenrich_Plots","_",main,"%03d.pdf",sep=""),onefile = F, width=14,height=14,pointsize = 8)
    #dev.off()
  } else {
    stop("Please specify a plot type !!")
  }
}

####################### Analysis started ###########################
# convert into character listofgene
listofgene = as.character(gene.target$mRNA)
# ############### conversion no needed if user already have entrez id!! ####################
# ############### no needed of conversion of id if user already have entrez id!! ####################
if (entrezID=="NO" || entrezID=="no" || entrezID=="No") {
  cat("#ID conversion started..\n")
  id <- suppressMessages(suppressWarnings(bitr(listofgene, fromType=GeneID, toType="ENTREZID", OrgDb=OrgDb, drop = TRUE)))
  ### convert as charecter vector
  id2 = as.character(id$ENTREZID)
  ## to strip the NA values
  x <- id2[!is.na(id2)] ## entrez id
  x <- unique(x)
  write.table(id, paste(output.dir.path,"/",sps,"_Converted_Gene_ID.txt",sep=""),quote=F,row.names = F)
  cat("#ID conversion Done.\n")
} else if (entrezID=="YES" || entrezID=="yes" || entrezID=="Yes") {
  ## No conversion needed
  cat ("#No conversion Needed.\n")
  x <- as.character(listofgene)
} else {
  stop("#Please check your input for: entrezID(yes/no), OrgDb and GeneID !!)")
}
###############################################################
######Calling function from {clusterProfiler}############

### Ontology should be DO if sps is hsa
if (sps=="hsa") {
  if (Ontology=="DO" || Ontology=="Do" || Ontology=="do") {
    cat("#Disease_Ontology analysis being processed..\n")
    enrich_go_do <- enrichDO(gene=as.vector(x), ont="DO", pvalueCutoff=as.numeric(pvalueCutoff),qvalueCutoff=as.numeric(qvalueCutoff), readable=TRUE, minGSSize = 5)
    ## print results
    print.go.out(goResults = enrich_go_do,outName = "_Disease_Ontology")
    cat("#Done.\n")
  }
} else if (sps!="hsa" && Ontology=="DO" || Ontology=="Do" || Ontology=="do") {
  stop("!! Disease ontology(DO) is only valid for Human (hsa); please check input for : sps and organism.\n")
}

##head(as.data.frame(enrich_go_do))
### gene ontology if ANY (ie.bp,cc and mf)
#keytype = c("ENTREZID", "ENSEMBL", "SYMBOL")
if (Ontology == "ANY" || Ontology == "any" || Ontology == "Any") {
  ### GO Enrichment analysis  from  {clusterProfiler}
  cat("#Performing Gene Ontology analysis...\n")
  message("Enrichment of Biological Process..")
  enrich_go_bp <- enrichGO(gene=as.vector(x), OrgDb=OrgDb, ont="BP", pvalueCutoff=as.numeric(pvalueCutoff), qvalueCutoff=as.numeric(qvalueCutoff), readable=TRUE, minGSSize = 5)
  message("Done.\n")
  message("Enrichment of Cellular Components..")
  enrich_go_cc <- enrichGO(gene=as.vector(x), OrgDb=OrgDb, ont="CC", pvalueCutoff=as.numeric(pvalueCutoff), qvalueCutoff=as.numeric(qvalueCutoff), readable=TRUE, minGSSize = 5)
  message("Done.\n")
  message("Enrichment of Molecular Functions..")
  enrich_go_mf <- enrichGO(gene=as.vector(x), OrgDb=OrgDb, ont="MF", pvalueCutoff=as.numeric(pvalueCutoff), qvalueCutoff=as.numeric(qvalueCutoff), readable=TRUE, minGSSize = 5)
  message("Done.")
  cat("#Analysis of GO successfully completed.\n")
}

# # summary
if (Ontology == "ANY" || Ontology == "any" || Ontology == "Any") {
  ## create GO enriched output
  cat("#Printing output of GO (Biological_Process) analysis..\n")
  print.go.out(goResults = enrich_go_bp, outName = "_Biological_Process")
  cat("Done.\n")
  cat("#Printing output of GO (Molecular_Function) analysis..\n")
  print.go.out(goResults = enrich_go_mf, outName = "_Molecular_Function")
  cat("Done.\n")
  cat("#Printing output of GO (Cellular_Compenents) analysis..\n")
  print.go.out(goResults = enrich_go_cc, outName = "_Cellular_Compenents")
  cat("Done.\n")
}

##define a null object for pathways
enrich_path <- NULL
## pathway analysis by Reactome from {ReactomePA}
if (pathways == "reactome" || pathways == "REACTOME" || pathways == "Reactome") {
  cat("#Predicting reactome pathways..\n")
  REACTOMEResult <- enrichPathway(gene=as.vector(x), organism = organism, pvalueCutoff = as.numeric(pvalueCutoff), pAdjustMethod = "BH", qvalueCutoff = as.numeric(qvalueCutoff), minGSSize = 5, readable=T)
  enrich_path <- REACTOMEResult
  cat("#Printing output REACTOMEResult...\n")
  print.go.out(goResults = REACTOMEResult, outName = "_ReactomePathways")
  cat("Done.\n")
}

### Check if internet accessible
if (download_latest_kegg == "no" || download_latest_kegg == "No" || download_latest_kegg == "NO") {
  download_latest = TRUE
} else {
  download_latest = FALSE ## in order to download new updated path from kegg website
}

## kegg pathways analysis
if (pathways=="kegg" || pathways=="Kegg" || pathways=="KEGG") {
  cat("#KEGG Enrichment Analysis is processed..\n")
  KEGGresult <- enrichKEGG(gene=as.vector(x), organism = as.character(sps), pvalueCutoff = as.numeric(pvalueCutoff), qvalueCutoff = as.numeric(qvalueCutoff), minGSSize = 5, use_internal_data = download_latest)
  enrich_path <- KEGGresult
  cat("#Printing output KEGGresult...\n")
  print.go.out(goResults = KEGGresult, outName = "_KEGG_Pathways")
  cat("done.\n")
}

## kegg pathways can be plotted but should have atleast one pathways information
## kegg pathway id can be view by:
# # row.names(as.data.frame(KEGGresult))
# # browseKEGG(KEGGresult, "hsa04142")

##=== function to plot enrichment
plotEnrichment <- function (enrichResults, Title) {
  ## call plot function type to plots network and GO enrichments
  print(barplot(enrichResults, drop=T, showCategory = 30, title = Title))
  # altrneative of barplot
  print(DOSE::dotplot(enrichResults, x = "geneRatio", colorBy = "p.adjust", showCategory = 10, title = Title))
  tl <- paste0(enrichResults@ontology,"_EnrichMap")
  print(DOSE::enrichMap(enrichResults,n = 10, main=tl))
}


## pathway analysis plots
if (pathways=="reactome" || pathways=="Reactome" || pathways=="REACTOME" || pathways=="kegg" || pathways=="Kegg" || pathways=="KEGG") {
  ## check enrich_path
  if (is.null(enrich_path)) {
    message("#Output of pathways analysis is empty!")
    message("#Try with different pathways !!")
    cat("#Tip: Please also make sure that you have internet connectivity in case of download_latest_kegg=yes.\n")
  } else if (nrow(as.data.frame(enrich_path)) >=5) {
    cairo_pdf(paste(output.dir.path,"/",sps,"_PathwaysEnrichment_Plots.pdf",sep=""),width=14,height=14,pointsize = 6)
    plotEnrichment (enrich_path, Title = "Ontology")
    graphics.off()
  } else {
    cat("#Number of predicted pathways are less than five or NULL, plottings skipping.\n")
    cat("#Tip: Use lower cut off of P-value and q-values to get some.\n")
    cat("#Tip: Please also make sure that you have internet connectivity in case of download_latest_kegg=yes.\n")
  }
}


## visual analysis of ontology (mf)
if (Ontology == "ANY" || Ontology == "any" || Ontology == "Any") {
  ## plot open
  if (nrow(as.data.frame(enrich_go_mf)) >=5) {
    plotingtype(plotType = plotType, main = "MF")
    message("#Plotting output of Molecular Function...\n")
    plotEnrichment (enrich_go_mf, Title="Molecular Function")
    message("done.\n")
    graphics.off()
  } else {
    message("#can not plot Molecular Function ontology, due to its null or less than five...\n")
  }
}

## visual analysis of ontology (cc)
if (Ontology == "ANY" || Ontology == "any" || Ontology == "Any") {
  ## plot open
  if (nrow(as.data.frame(enrich_go_cc)) >=5) {
    plotingtype(plotType = plotType, main = "CC")
    message("#Plotting output of Cellular Component...\n")
    plotEnrichment (enrich_go_cc, Title="Cellular Component")
    message("done.\n")
    graphics.off()
  } else {
    message("#can not plot Cellular Component ontology, due to its null or less than five...\n")
  }
}

## visual analysis of ontology (bp)
if (Ontology == "ANY" || Ontology == "any" || Ontology == "Any") {
  ## plot open
  if (nrow(as.data.frame(enrich_go_bp)) >=5) {
    plotingtype(plotType = plotType, main="BP")
    message("#Plotting output of Biological Process...\n")
    plotEnrichment (enrich_go_bp, Title="Biological Process")
    message("done.\n")
    graphics.off()
  } else {
    message("#can not plot Biological Process ontology, due to its null or less than five...\n")
  }
}


## visual analysis of ontology (do)
if (Ontology=="DO" || Ontology=="Do" || Ontology=="do") {
  ## plot open
  if (nrow(as.data.frame(enrich_go_do)) >=5) {
    plotingtype(plotType = plotType, main="DO")
    message("#Plotting output of Disease Ontology...\n")
    plotEnrichment (enrich_go_do, Title="Disease Ontology")
    message("done.\n")
    graphics.off()
  } else {
    message("#can not plot Disease Ontologyy, due to its null or less than five...\n")
  }
}


##dev.off()
#########################################################
## save work space
cat ("#Saving Rdata....\n")
save.image(paste(output.dir.path,"/",sps,"_GeneOntologyNetwork.RData", sep=""))
cat("done.\n")
cat("\n#Gene Ontology and Network analysis has successfully performed, go to your output.dir.path to check results.\n\n")
