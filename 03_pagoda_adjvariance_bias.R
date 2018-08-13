#!/bin/Rscript

###########################
### REQUIRED PACKAGES
###########################

library(optparse)
library(scde)
library(biomaRt)
library(ShortRead)
library(RColorBrewer)
library(coin)
library(ggplot2)
library(reshape2)
library(plyr)
library(GGally)

###########################
### FUNCTIONS
###########################

###########################
### GLOBALS
###########################

#Main paths
base_path<-getwd() #Assumes current working directory is repository base folder ("chromatin_noise_paper")
misc_functions_path<-file.path(base_path, 'misc_functions.R')

option_list <- list(
  make_option(opt_str=c("--inputFile", "-i"), help = "Path to pagoda_varnorm RData file containing normalised count data"),
  make_option(opt_str=c("--geneSetFile", "-e"), help = "Path to RData file containing gene set data"),
  make_option(opt_str=c("--GOFile", "-g"), help = "Path to Gene Ontology data"),
  make_option(opt_str=c("--numCores", "-c"), type='integer', help = "Number of CPU cores"),
  make_option(opt_str=c("--nRandomizations", "-n"), type='integer', help = "Number of random gene sets (of the same size) to be evaluated in parallel with each gene set"),
  make_option(opt_str=c("--numComponents", "-t"), type='integer', default=1, help = "Number of principal components to determine for each gene set"),
  make_option(opt_str=c("--maxPathwaySize", "-m"), type='integer', help = "Maximum number of observed genes in a valid gene set"),
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

#Load gene set data
load(arg_list$geneSetFile)

#Load GO data
load(arg_list$GOFile)

#Required functions
source(misc_functions_path)

###########################

#Remove irrelevant GO/genes and transform to environment
for(i in names(GOBPgene_3ontologies)){
  GOBPgene_3ontologies[[i]]<-GOBPgene_3ontologies[[i]][GOBPgene_3ontologies[[i]] %in% names(genevar)]
}
GOBPgene_3ontologies<-GOBPgene_3ontologies[sapply(GOBPgene_3ontologies, length)!=0]
go.env_GOBP<-list2env(GOBPgene_3ontologies)
for(i in names(GOCCgene_3ontologies)){
  GOCCgene_3ontologies[[i]]<-GOCCgene_3ontologies[[i]][GOCCgene_3ontologies[[i]] %in% names(genevar)]
}
GOCCgene_3ontologies<-GOCCgene_3ontologies[sapply(GOCCgene_3ontologies, length)!=0]
go.env_GOCC<-list2env(GOCCgene_3ontologies)
for(i in names(GOMFgene_3ontologies)){
  GOMFgene_3ontologies[[i]]<-GOMFgene_3ontologies[[i]][GOMFgene_3ontologies[[i]] %in% names(genevar)]
}
GOMFgene_3ontologies<-GOMFgene_3ontologies[sapply(GOMFgene_3ontologies, length)!=0]
go.env_GOMF<-list2env(GOMFgene_3ontologies)

###########################

#Mask out genes with length <10kb for whole gene interval overlap features
temp_excluded<-colnames(elementMetadata(allgeneInfo_gr))[c(grep("globalchromatin.* band", colnames(elementMetadata(allgeneInfo_gr))),
                                                           grep("globalchromatin.*-chromosome", colnames(elementMetadata(allgeneInfo_gr))),
                                                           grep("globalchromatin_Super", colnames(elementMetadata(allgeneInfo_gr))),
                                                           grep("sequence_mRNA", colnames(elementMetadata(allgeneInfo_gr))),
                                                           grep("sequence_NMD repressed", colnames(elementMetadata(allgeneInfo_gr))),
                                                           grep("sequence_Short", colnames(elementMetadata(allgeneInfo_gr))),
                                                           grep("localchromatin.* target", colnames(elementMetadata(allgeneInfo_gr))))]
allgeneInfo_gr<-mask_minlength(allgeneInfo_gr, min_length=1e4, feature_names_excluded=temp_excluded)

###########################

#Rank-based gene expression variation measure
genevar_rank<-binx_nearesty(varinfo$avmodes, genevar[names(varinfo$avmodes)], bin_size=101)
#Remove NAs
genevar_rank<-genevar_rank[!is.na(genevar_rank)]

###########################

#Test biased variance of all chromatin gene sets
gs_tests_all<-list()
gs_tests_all[['tssInfo']]<-genesets_biased(as.list(elementMetadata(tssInfo_gr[names(genevar_rank),])), genevar_rank)
gs_tests_all[['corePromInfo']]<-genesets_biased(as.list(elementMetadata(corePromInfo_gr[names(genevar_rank),])), genevar_rank)
gs_tests_all[['promInfo']]<-genesets_biased(as.list(elementMetadata(promInfo_gr[names(genevar_rank),])), genevar_rank)
gs_tests_all[['allgeneInfo']]<-genesets_biased(as.list(elementMetadata(allgeneInfo_gr[names(genevar_rank),])), genevar_rank)
gs_tests_all<-do.call('rbind', gs_tests_all)

###########################

#Test biased variance of all GO gene sets
gs_tests_all_GO<-list()
gs_tests_all_GO[['biological_process']]<-genesets_biased(GOBPgene_3ontologies, genevar_rank)
gs_tests_all_GO[['cellular_component']]<-genesets_biased(GOCCgene_3ontologies, genevar_rank)
gs_tests_all_GO[['molecular_function']]<-genesets_biased(GOMFgene_3ontologies, genevar_rank)

###########################

#Evaluate overdispersion of all gene sets
env_all<-list2env(c(binaryvect_to_geneset(allgeneInfo_gr[names(allgeneInfo_gr) %in% names(varinfo$avmodes),], 'allgeneInfo'),
                    binaryvect_to_geneset(promInfo_gr[names(promInfo_gr) %in% names(varinfo$avmodes),], 'promInfo'),
                    binaryvect_to_geneset(corePromInfo_gr[names(corePromInfo_gr) %in% names(varinfo$avmodes),], 'corePromInfo'),
                    binaryvect_to_geneset(tssInfo_gr[names(tssInfo_gr) %in% names(varinfo$avmodes),], 'tssInfo'),
                    GOBPgene_3ontologies, GOCCgene_3ontologies, GOMFgene_3ontologies))
#Calculate weighted first principal component magnitudes for each gene set
set.seed(1)
pwpca<-pagoda.pathway.wPCA(varinfo, env_all, n.components = arg_list$numComponents, n.cores = arg_list$numCores, n.randomizations=arg_list$nRandomizations, max.pathway.size = arg_list$maxPathwaySize)
pwpca_topdf<-pagoda.top.aspects(pwpca, return.table = TRUE, plot = F, z.score = -Inf, adjust.scores=F)

###########################

#Reformat results
rownames(pwpca_topdf)<-pwpca_topdf$name
pwpca_topdf$effect_size<-pwpca_topdf$score
pwpca_topdf$effect_size[grep('_no_', rownames(pwpca_topdf))]<-(-1)*pwpca_topdf$effect_size[grep('_no_', rownames(pwpca_topdf))]
rownames(pwpca_topdf)<-gsub('Info_yes_', 'Info.', rownames(pwpca_topdf))
rownames(pwpca_topdf)<-gsub('Info_no_', 'Info.', rownames(pwpca_topdf))
pwpca_topdf$p_value<-pnorm(pwpca_topdf$z, lower.tail=F)*2
pwpca_topdf$FDR<-p.adjust(pwpca_topdf$p_value, method='BH')

#Final chromatin gene set results
gsr_chr<-pwpca_topdf[grep('Info.', rownames(pwpca_topdf)),]
#Feauture category
gsr_chr$feature_category<-extract_feature_category(rownames(gsr_chr))
#Feauture names
gsr_chr$feature_name<-extract_feature_name(rownames(gsr_chr))
gsr_chr$set_covariability<-log2(abs(gsr_chr$effect_size))
#Biased variance effect size
gsr_chr$set_variability<-gs_tests_all[rownames(gsr_chr),]$effect_size
#Set sign of biased variance effect size accordingly
gsr_chr$set_variability[grep('_no_', gsr_chr$name)]<-1-gsr_chr$set_variability[grep('_no_', gsr_chr$name)]
#Gene region
gsr_chr$region<-as.factor(sapply(strsplit(rownames(gsr_chr), '\\.'), '[', 1))
#Gene region
gsr_chr$set_status<-'in'
gsr_chr$set_status[grep('_no_', gsr_chr$name)]<-'out'
gsr_chr$set_status<-as.factor(gsr_chr$set_status)

#Final GO gene set results
gsr_go<-pwpca_topdf[-grep('Info.', rownames(pwpca_topdf)),]
#Feature category
gsr_go$feature_category<-GOinfo_3ontologies[rownames(gsr_go),]$namespace_1003
#Feauture names
gsr_go$feature_name<-GOinfo_3ontologies[rownames(gsr_go),]$name_1006
gsr_go$set_covariability<-log2(abs(gsr_go$effect_size))
#Biased variance effect size
gsr_go$set_variability<-NA
gsr_go$set_variability[gsr_go$feature_category=='biological_process']<-gs_tests_all_GO[['biological_process']][rownames(gsr_go)[gsr_go$feature_category=='biological_process'],]$effect_size
gsr_go$set_variability[gsr_go$feature_category=='cellular_component']<-gs_tests_all_GO[['cellular_component']][rownames(gsr_go)[gsr_go$feature_category=='cellular_component'],]$effect_size
gsr_go$set_variability[gsr_go$feature_category=='molecular_function']<-gs_tests_all_GO[['molecular_function']][rownames(gsr_go)[gsr_go$feature_category=='molecular_function'],]$effect_size

###########################

#Save objects
save(countsMmus, countsERCC, cd, knn, varinfo, genevar, genevar_rank,
     allgeneInfo_gr, promInfo_gr, corePromInfo_gr, tssInfo_gr, 
     gs_tests_all, gs_tests_all_GO,
     env_all, pwpca, gsr_chr, gsr_go, file=paste(arg_list$projectName, '.RData', sep=''))

###########################
