#!/bin/Rscript

###########################
### REQUIRED PACKAGES
###########################

library(optparse)
library(scde)
library(biomaRt)
library(ShortRead)
library(ggplot2)
library(reshape2)
library(plyr)
library(GGally)

###########################
### FUNCTIONS
###########################

gene_pattern_heatmap<-function(input_pattern, output_file, width=10, height=4, colour_clip=4){
  plot_df<-melt(varinfo$mat[input_pattern$lab,order(input_pattern$oc)])
  colnames(plot_df)<-c('gene_id', 'cell', 'pattern')
  plot_df$pattern[plot_df$pattern>colour_clip]<-colour_clip
  plot_df$pattern[plot_df$pattern<(-colour_clip)]<-(-colour_clip)
  plot_df$gene_name<-factor(allgenes_info[as.character(plot_df$gene_id),1], levels=unique(allgenes_info[as.character(plot_df$gene_id),1]))
  p <- ggplot(plot_df, aes(cell, gene_name)) + geom_tile(aes(fill = pattern), colour = "white") + theme(axis.text.x = element_text(angle = 330, hjust = 0, size = 5))
  p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)
  ggsave(file=output_file, width=width, height=height)
}

###########################
### GLOBALS
###########################

option_list <- list(
  make_option(opt_str=c("--inputFile", "-i"), help = "Path to pagoda_varnorm RData file containing normalised count data"),
  make_option(opt_str=c("--GOFile", "-g"), help = "Path to Gene Ontology data"),
  make_option(opt_str=c("--numCores", "-c"), type='integer', help = "Number of CPU cores"),
  make_option(opt_str=c("--nRandomizations", "-n"), type='integer', help = "Number of random gene sets (of the same size) to be evaluated in parallel with each gene set"),
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

#Load GO data
load(arg_list$GOFile)

###########################

#Remove irrelevant GO/genes and transform to environment
for(i in names(GOBPgene_3ontologies)){
  GOBPgene_3ontologies[[i]]<-GOBPgene_3ontologies[[i]][GOBPgene_3ontologies[[i]] %in% names(genevar)]
}
GOBPgene_3ontologies<-GOBPgene_3ontologies[sapply(GOBPgene_3ontologies, length)!=0]
go.env_GOBP<-list2env(GOBPgene_3ontologies)

###########################

#Iterately remove cell cycle pattern until non-significant source of variation
correct_iter<-1
while(TRUE){
  print(correct_iter)
  #Calculate weighted first principal component magnitudes for each gene set (BP)
  set.seed(1)
  pwpca_GOBP<-pagoda.pathway.wPCA(varinfo, go.env_GOBP, n.components = 1, n.cores = arg_list$numCores, n.randomizations=arg_list$nRandomizations, max.pathway.size = arg_list$maxPathwaySize)
  #Get all aspects
  df_GOBP<-pagoda.top.aspects(pwpca_GOBP, return.table = TRUE, z.score = -Inf, adjust.scores=T)
  #Plot cell cycle pattern
  temp_name<-rownames(GOinfo_3ontologies)[GOinfo_3ontologies$name_1006=='cell cycle']
  temp_patterndetails<-pagoda.show.pathways(temp_name, varinfo, go.env_GOBP, showRowLabels = TRUE, return.details=T, n.genes=20, plot=F)
  gene_pattern_heatmap(temp_patterndetails, paste('pagoda_cellcycle_pattern_allgenes_heatmap_iteration', correct_iter, '.pdf', sep=''), width=10, height=5)
  #Significance of cell cycle pattern
  temp_sig<-pnorm(df_GOBP[df_GOBP$name==temp_name,]$adj.z, lower.tail=F)*2
  if(temp_sig<0.05){
    print(temp_sig)
    #Subtract the pattern
    varinfo<-pagoda.subtract.aspect(varinfo, temp_patterndetails$scores[,1])      
  }else{
    break
  }
  correct_iter<-correct_iter+1
}

#Recalculate adjusted variance
genevar <- scde:::weightedMatVar(varinfo$mat, varinfo$matw)

###########################

#Save objects
save(countsMmus, countsERCC, cd, knn, varinfo, genevar, file=paste(arg_list$projectName, '.RData', sep=''))

###########################


