#!/bin/Rscript

###########################
### REQUIRED PACKAGES
###########################

library(optparse)
library(scde)
library(biomaRt)
library(ShortRead)

###########################
### FUNCTIONS
###########################

###########################
### GLOBALS
###########################

option_list <- list(
  make_option(opt_str=c("--inputFile", "-i"), help = "Path to preQC RData file containing count data and (optionally) batch description"),
  make_option(opt_str=c("--maxTranscriptLengthFile", "-x"), help = "Path to file with length of longest possible transcript (transcript_maxlength)"),
  make_option(opt_str=c("--spikeIn", "-k"), type='logical', default=TRUE, help = "ERCC spike-ins present"),
  make_option(opt_str=c("--numCores", "-c"), type='integer', help = "Number of CPU cores"),
  make_option(opt_str=c("--minDetected", "-d"), type='integer', help = "Minimum number of cells a gene must be seen in. Genes not seen in a sufficient number of cells will be removed"),
  make_option(opt_str=c("--minCountThreshold", "-m"), type='integer', help = "Minimum number of reads required for a measurement to be considered non-failed"),
  make_option(opt_str=c("--maxAdjVar", "-a"), type='integer', help = "maximum value allowed for the estimated adjusted variance (capping of adjusted variance is recommended when scoring pathway overdispersion relative to randomly sampled gene sets)"),
  make_option(opt_str=c("--biomartSpecies", "-s"), default='mmusculus', help = "BiomaRt species names"),
  make_option(opt_str=c("--biomartHost", "-t"), default='apr2013.archive.ensembl.org', help = "BiomaRt host URL"),
  make_option(opt_str=c("--chromosomeSubset", "-r"), default=paste(c(1:19, 'X', 'Y', 'MT'), collapse=','), help = "Comma-separated list of (Ensembl) chromosomes of interest"),
  make_option(opt_str=c("--outputPath", "-o"), help = "Path to directory to use for output files"),
  make_option(opt_str=c("--projectName", "-p"), help = "Project name")
)

arg_list<-parse_args(OptionParser(option_list=option_list))

###########################
### MAIN
###########################

setwd(arg_list$outputPath)

#Load input data
load(arg_list$inputFile)

###########################

#Separate count data
countsMmus<-NULL
countsERCC<-NULL
if(arg_list$spikeIn){
  countsMmus<-dataSC_filter[-which(substr(rownames(dataSC_filter), 1, 4)=="ERCC"),]
  countsERCC<-dataSC_filter[which(substr(rownames(dataSC_filter), 1, 4)=="ERCC"),]    
}else{
  countsMmus<-dataSC_filter
}

#Remove poor cells and genes
cd<-clean.counts(countsMmus, min.detected=arg_list$minDetected)

#Fitting error models
knn<-knn.error.models(cd, k = as.integer(ncol(cd)/4), n.cores = arg_list$numCores, min.count.threshold = arg_list$minCountThreshold, max.model.plots = 10)

#Control for gene length if not UMI-based dataset
load(arg_list$maxTranscriptLengthFile)
all_lengths<-transcript_maxlength[rownames(cd)]
names(all_lengths)<-rownames(cd)

#Control for batch if present
all_batch<-NULL
if('dataSC_filter_batch' %in% ls()){
  names(dataSC_filter_batch)<-colnames(dataSC_filter)
  all_batch<-dataSC_filter_batch[colnames(cd)]
}

#Normalizing variance
pdf('normalizing_variance.pdf', width=10, height=5)
varinfo<-pagoda.varnorm(knn, counts = cd, max.adj.var = arg_list$maxAdjVar, n.cores = arg_list$numCores, plot = TRUE, gene.length=all_lengths, batch=all_batch)
dev.off()

#Controlling for sequencing depth
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))
#Recalculate adjusted variance
genevar <- scde:::weightedMatVar(varinfo$mat, varinfo$matw)

###########################

#Save objects
save(countsMmus, countsERCC, cd, knn, varinfo, genevar, file=paste(arg_list$projectName, '.RData', sep=''))

###########################
