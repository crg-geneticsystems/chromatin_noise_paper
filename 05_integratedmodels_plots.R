
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
library(DESeq2)
library(rtracklayer)
library(randomForest)
library(cluster)
library(ppcor)
library(plyr)
library(glmnet)

###########################
### FUNCTIONS
###########################

#Integrated model wrapper
model_wrapper<-function(dataset_name, model_type='glm', feature_subset='coreprom', exclude_mid_tercile=T){
  temp<-NULL
  if(feature_subset=='coreprom'){
    #Only include core promoter features
    temp<-as.data.frame(elementMetadata(corePromInfo_gr))
    rownames(temp)<-names(corePromInfo_gr)
    #Subset features
    temp<-temp[,1:10]
    temp<-temp[,-4]
    temp<-temp[names(genevarrank_list[[dataset_name]]),]
  }else if(feature_subset=='coreprom_sharpbroad'){
    #Only include core promoter features
    temp<-as.data.frame(elementMetadata(corePromInfo_gr))
    rownames(temp)<-names(corePromInfo_gr)
    #Subset features
    temp<-temp[,1:10]
    # #Exclude ambiguous sharp/broad regions
    # temp$sequence_Sharp.TSS[which(temp$sequence_Sharp.TSS & temp$sequence_Broad.TSS)]<-NA
    temp<-temp[,-c(4, 8)]
    temp<-temp[names(genevarrank_list[[dataset_name]]),]
  }else if(feature_subset=='coreprom_promwidth'){
    #Only include core promoter features
    temp<-as.data.frame(elementMetadata(corePromInfo_gr))
    rownames(temp)<-names(corePromInfo_gr)
    #Subset features
    temp<-temp[,1:10]
    temp<-temp[,-c(4, 8, 9)]
    temp$sequence_Prom.width<-log10(TSSclust_width[rownames(temp)])
    temp<-temp[names(genevarrank_list[[dataset_name]]),]
  }else if(feature_subset=='coreprom_chrom'){
    #Include core promoter features, chromatin features (include LADs)
    temp<-temp_nnfreq[names(genevarrank_list[[dataset_name]]),]
    temp$sequence_CpG_island<-corePromInfo_gr[rownames(temp),]$"sequence_CpG island"
    temp$sequence_TATAbox<-as.numeric(coreProm_TATAscores[rownames(temp)]>FPR_10)
    temp$allgeneInfo_cLAD<-as.numeric(countOverlaps(allgeneInfo_gr[rownames(temp),], lad_gr)!=0)
  }else if(feature_subset=='coreprom_chrom_cLADNA'){
    #Include core promoter features, chromatin features (include LADs)
    temp<-temp_nnfreq[names(genevarrank_list[[dataset_name]]),]
    temp$sequence_CpG_island<-corePromInfo_gr[rownames(temp),]$"sequence_CpG island"
    temp$sequence_TATAbox<-as.numeric(coreProm_TATAscores[rownames(temp)]>FPR_10)
    temp$allgeneInfo_cLAD<-allgeneInfo_gr[rownames(temp),]$"globalchromatin_cLAD"
  }else if(feature_subset=='coreprom_chrom_X'){
    #Include core promoter features, chromatin features (include LADs and chrX)
    temp<-temp_nnfreq[names(genevarrank_list[[dataset_name]]),]
    temp$sequence_CpG_island<-corePromInfo_gr[rownames(temp),]$"sequence_CpG island"
    temp$sequence_TATAbox<-as.numeric(coreProm_TATAscores[rownames(temp)]>FPR_10)
    temp$allgeneInfo_cLAD<-as.numeric(countOverlaps(allgeneInfo_gr[rownames(temp),], lad_gr)!=0)
    temp$allgeneInfo_chrX<-allgeneInfo_gr[rownames(temp),]$"globalchromatin_X-chromosome"
  }else if(feature_subset=='coreprom_chrom_X_mRNAstability'){
    #Include core promoter features, chromatin features (include LADs and chrX and mRNA stability)
    temp<-temp_nnfreq[names(genevarrank_list[[dataset_name]]),]
    temp$sequence_CpG_island<-corePromInfo_gr[rownames(temp),]$"sequence_CpG island"
    temp$sequence_TATAbox<-as.numeric(coreProm_TATAscores[rownames(temp)]>FPR_10)
    temp$allgeneInfo_cLAD<-as.numeric(countOverlaps(allgeneInfo_gr[rownames(temp),], lad_gr)!=0)
    temp$allgeneInfo_chrX<-allgeneInfo_gr[rownames(temp),]$"globalchromatin_X-chromosome"
    temp$mRNA_stability<-mRNA_decay[rownames(temp)]
  }else if(feature_subset=='coreprom_chrom_X_mRNAstability_SE'){
    #Include core promoter features, chromatin features (include LADs and chrX and mRNA stability and SE)
    temp<-temp_nnfreq[names(genevarrank_list[[dataset_name]]),]
    temp$sequence_CpG_island<-corePromInfo_gr[rownames(temp),]$"sequence_CpG island"
    temp$sequence_TATAbox<-as.numeric(coreProm_TATAscores[rownames(temp)]>FPR_10)
    temp$allgeneInfo_cLAD<-as.numeric(countOverlaps(allgeneInfo_gr[rownames(temp),], lad_gr)!=0)
    temp$allgeneInfo_chrX<-allgeneInfo_gr[rownames(temp),]$"globalchromatin_X-chromosome"
    temp$allgeneInfo_SE<-allgeneInfo_gr[rownames(temp),]$"globalchromatin_Super-enhancer ES (target)"
    temp$mRNA_stability<-mRNA_decay[rownames(temp)]
  }else if(feature_subset=='coreprom_chrom_restrict'){
    #Include core promoter features, chromatin features (include LADs)
    temp<-temp_nnfreq[names(genevarrank_list[[dataset_name]]),]
    temp$sequence_CpG_island<-corePromInfo_gr[rownames(temp),]$"sequence_CpG island"
    temp$sequence_TATAbox<-corePromInfo_gr[rownames(temp),]$"sequence_TATAbox"
    temp$allgeneInfo_cLAD<-as.numeric(countOverlaps(allgeneInfo_gr[rownames(temp),], lad_gr)!=0)
  }else if(feature_subset=='coreprom_chrom_X_mRNAstability_SE_restrict'){
    #Include core promoter features, chromatin features (include LADs and chrX and mRNA stability and SE)
    temp<-temp_nnfreq[names(genevarrank_list[[dataset_name]]),]
    temp$sequence_CpG_island<-corePromInfo_gr[rownames(temp),]$"sequence_CpG island"
    temp$sequence_TATAbox<-corePromInfo_gr[rownames(temp),]$"sequence_TATAbox"
    temp$allgeneInfo_cLAD<-as.numeric(countOverlaps(allgeneInfo_gr[rownames(temp),], lad_gr)!=0)
    temp$allgeneInfo_chrX<-allgeneInfo_gr[rownames(temp),]$"globalchromatin_X-chromosome"
    temp$allgeneInfo_SE<-allgeneInfo_gr[rownames(temp),]$"globalchromatin_Super-enhancer ES (target)"
    temp$mRNA_stability<-mRNA_decay[rownames(temp)]
  }
  temp$expr_noise<-genevarrank_list[[dataset_name]]
  temp$gene_length<-log10(width(allgeneInfo_gr[names(genevarrank_list[[dataset_name]]),]))
  temp$mean_expr<-log10(meansMmus_list[[dataset_name]][names(genevarrank_list[[dataset_name]])])
  #Subset rows
  temp<-temp[!is.na(apply(temp, 1, sum)),]
  print(dim(temp)[1])
  if(exclude_mid_tercile){
    #Remove middle tercile
    temp<-temp[temp$expr_noise<=34 | temp$expr_noise>=68,]
  }
  # if(exclude_mid_tercile){
  #   #Remove middle tercile
  #   temp_tercile<-quantile(temp$expr_noise, c(1/3, 2/3))
  #   temp_tercile_num<-sum(temp$expr_noise<=temp_tercile[1] | temp$expr_noise>temp_tercile[2])
  #   temp<-temp[temp$expr_noise<=temp_tercile[1] | temp$expr_noise>temp_tercile[2],]
  # }
  #Remove random tercile
  #   temp<-temp[sample(1:dim(temp)[1], temp_tercile_num),]
  #Subset columns
  #   temp<-temp[,-which(apply(temp, 2, sum)<20)]
  #Discretise noise measure
  temp$expr_noise<-as.factor(as.numeric(temp$expr_noise>51))
  #Normalise
  temp<-glm_prenorm(temp)
  #Model
  int_model<-as.data.frame(summary(glm(expr_noise~., family=binomial(link='logit'), data=temp))$coefficients)[-1,]
  if(model_type=='glm'){
    int_model$dataset = dataset_name
    #glm
    return(int_model)
  }else if(model_type=='glm_uni'){
    #Univariate models using glm
    uni_model_list<-list()
    for(i in 1:dim(temp)[2]){
      temp_name<-colnames(temp)[i]
      if(temp_name!='expr_noise'){
        uni_model_list[[i]]<-summary(glm(expr_noise~., family=binomial(link='logit'), data=temp[,colnames(temp) %in% c('expr_noise', temp_name)]))
      }
    }
    uni_model<-do.call('rbind', sapply(sapply(uni_model_list, '[', 'coefficients'), as.data.frame))
    uni_model<-uni_model[-grep('Intercept', rownames(uni_model)),]
    colnames(uni_model)<-colnames(int_model)
    uni_model$dataset = dataset_name
    return(uni_model)
  }else if(model_type=='glmnet'){
    #glmnet
    temp_glmnet<-sapply(temp, unlist)
    temp_glmnet_out<-cv.glmnet(x=temp_glmnet[,-grep('expr_noise', colnames(temp_glmnet))], y=temp_glmnet[,"expr_noise"], standardize=F, family="binomial")
    glmnet_model<-as.data.frame(as.matrix(coef(temp_glmnet_out, s="lambda.1se")))
    colnames(glmnet_model)[1]<-'Estimate'
    glmnet_model$'Pr(>|z|)'<-as.numeric(glmnet_model$Estimate==0)
    glmnet_model$dataset = dataset_name
    return(glmnet_model[-1,])
  }
}

#Integrated model coefficient plotter
plot_coefficients<-function(model_1, model_2, model_3, variable_levels=NULL, output_file, width=5, height=5){
  #Plot coefficients
  plot_df_1<-data.frame(coefficient=model_1[,'Estimate'], pval=model_1[,'Pr(>|z|)'], feature=gsub('Info', '', gsub('_' , ' ', gsub('\\.', ' ', gsub('sequence_', '', rownames(model_1))))), dataset = model_1[,'dataset'], stringsAsFactors=F, row.names = rownames(model_1))
  plot_df_2<-data.frame(coefficient=model_2[,'Estimate'], pval=model_2[,'Pr(>|z|)'], feature=gsub('Info', '', gsub('_' , ' ', gsub('\\.', ' ', gsub('sequence_', '', rownames(model_2))))), dataset = model_2[,'dataset'], stringsAsFactors=F, row.names = rownames(model_2))
  plot_df_3<-data.frame(coefficient=model_3[,'Estimate'], pval=model_3[,'Pr(>|z|)'], feature=gsub('Info', '', gsub('_' , ' ', gsub('\\.', ' ', gsub('sequence_', '', rownames(model_3))))), dataset = model_3[,'dataset'], stringsAsFactors=F, row.names = rownames(model_3))
  plot_df<-rbind(plot_df_1, plot_df_2, plot_df_3)
  #Set non-signficant gene sets to open dot
  plot_df$significance<-'P<0.05'
  plot_df$significance[plot_df$pval>=0.05]<-'P>=0.05'
  plot_df$significance<-paste(plot_df$significance, plot_df$dataset, sep='.')
  plot_df$significance<-as.factor(plot_df$significance)
  plot_shapes<-c(1, 16, 2, 17, 0, 15, 16, 17, 15)
  names(plot_shapes)<-c('P>=0.05.islam','P<0.05.islam',
                        'P>=0.05.grunS', 'P<0.05.grunS',
                        'P>=0.05.grun2i', 'P<0.05.grun2i', 
                        'islam', 'grunS', 'grun2i')
  #Bubble colour
  plot_df$location<-'5 other'
  plot_df$location[grep('tssInfo_', rownames(plot_df))]<-'1 TSS'
  plot_df$location[grep('sequence_', rownames(plot_df))]<-'2 core promoter'
  plot_df$location[grep('corePromInfo', rownames(plot_df))]<-'2 core promoter'
  plot_df$location[grep('promInfo', rownames(plot_df))]<-'3 promoter'
  plot_df$location[grep('allgeneInfo', rownames(plot_df))]<-'4 whole gene body'
  plot_df$location<-as.character(plot_df$location)
  plot_cols<-c(brewer.pal(5, 'Spectral')[-3], 'grey')
  names(plot_cols)<-c('1 TSS', '2 core promoter', '3 promoter', '4 whole gene body', '5 other')
  #Plot order
  if(is.null(variable_levels)){
    mean_coeff<-ddply(plot_df, c("feature"), summarise, dfreq = mean(coefficient))
    sum_sig<-ddply(plot_df, c("feature"), summarise, dfreq = sum(pval<0.05))
    variable_levels<-as.character(mean_coeff$feature)[order(sign(mean_coeff$dfreq)*sum_sig$dfreq, mean_coeff$dfreq, decreasing=T)]
  }
  plot_df$feature<-factor(plot_df$feature, levels=variable_levels)
  if(dim(model_1)[2]==5){
    #glm or glm_uni results
    d <- ggplot(plot_df, aes(feature, coefficient, group = feature)) +
      #   geom_boxplot() + 
      geom_line(aes(group = feature, colour=location), position = position_dodge(width = 0.75)) +
      geom_point(aes(shape=dataset), colour='white', size=2, position = position_dodge(width = 0.75)) +
      geom_point(data=plot_df[as.character(plot_df$significance)=='P<0.05',], aes(colour=location, shape=significance), size=2, position = position_dodge(width = 0.75)) +
      geom_point(data=plot_df[as.character(plot_df$significance)!='P<0.05',], aes(colour=location, shape=significance), size=2, position = position_dodge(width = 0.75)) +
      scale_colour_manual(values = plot_cols) + 
      scale_fill_manual(values = plot_cols) + 
      scale_shape_manual(values = plot_shapes) + 
      theme_bw() + 
      geom_hline(yintercept = 0, size=0.5) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
      # coord_cartesian(ylim = c(-0.75, 1.5))
    ggsave(file=output_file, width=width, height=height)
    return(variable_levels)
  }else{
    #glmnet results
    d <- ggplot(plot_df, aes(feature, coefficient, group = feature)) +
      #   geom_boxplot() + 
      geom_line(aes(group = feature, colour=location), position = position_dodge(width = 0.75)) +
      # geom_point(color='white', size=2, position = position_dodge(width = 0.75)) +
      # geom_point(data=plot_df[as.character(plot_df$significance)=='P<0.05',], aes(colour=location, shape=significance), size=2, position = position_dodge(width = 0.75)) +
      # geom_point(data=plot_df[as.character(plot_df$significance)!='P<0.05',], aes(colour=location, shape=significance), size=2, position = position_dodge(width = 0.75)) +
      geom_point(aes(colour=location, shape=dataset), size=2, position = position_dodge(width = 0.75)) +
      scale_colour_manual(values = plot_cols) + 
      scale_fill_manual(values = plot_cols) + 
      # scale_shape_manual(values = plot_shapes) + 
      theme_bw() + 
      geom_hline(yintercept = 0, size=0.5) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
      # coord_cartesian(ylim = c(-0.75, 1.5))
    ggsave(file=output_file, width=width, height=height)
    return(variable_levels)
  }
}

###########################
### GLOBALS
###########################

#Main paths
base_path<-getwd() #Assumes current working directory is repository base folder ("chromatin_noise_paper")
misc_functions_path<-file.path(base_path, 'misc_functions.R')
misc_data_path<-file.path(base_path, 'misc_data_files')
project_name<-'05_integratedmodels_plots'
working_directory<-file.path(base_path, project_name)

#Gene annotations and gene set files
ensapi_file<-file.path(misc_data_path, 'ens71_allgenefeatures.txt')

#Miscellaneous data files
chip_sub_folder<-file.path(misc_data_path, "bam_bed_enrich_nobruce4input")
SEcor_file<-file.path(misc_data_path, 'cor_list_SE.RData')
TATAscore_file<-file.path(misc_data_path, 'SC_TATAscore.RData')
input_file_cLAD<-file.path(misc_data_path, 'GSE17051_HMM_state_calls_per_probe_contiguous_cLAD_mm10_UCSCliftover.bed')
mrna_decay_path_friedel<-file.path(misc_data_path, 'Friedel2009_NIH3T3_mRNAhalflife_unix.txt')
cage_consensus_clusters_path<-file.path(misc_data_path, 'CAGE_wbembryo_consensus_clusters_mm10_UCSCliftover.bed')

#Plot colour palette
palette_rdbu<-brewer.pal(7, 'RdBu')
col_red<-palette_rdbu[1]
col_blue<-palette_rdbu[7]
palette_brgr<-brewer.pal(7, 'BrBG')
col_brown<-palette_brgr[1]
col_green<-palette_brgr[7]

###########################
### MAIN
###########################

#Create output directory
dir.create(working_directory, showWarnings = FALSE)
setwd(working_directory)

#Required functions
source(misc_functions_path)

###########################
### Load required data and annotations
###########################

#Load input data (after cell cycle correction)
ncd_list<-list()
nCountsMmus_list<-list()
nCountsMmus_mean_list<-list()
meansMmus_list<-list()
genevarrank_list<-list()
temp_dirs<-list.dirs('..', recursive=F)
temp_dirs<-temp_dirs[grep('pagoda_adjvariance_bias', temp_dirs)]
#Retain cell cycle corrected data
temp_dirs<-temp_dirs[grep('ccc', temp_dirs)]
# temp_dirs<-temp_dirs[-grep("covariancecontrol|graf|klein|zeisel", temp_dirs)]
#Exclude pool-and-split samples
temp_dirs<-temp_dirs[-grep('PandS', temp_dirs)]
#Selected UMI-based datasets
temp_dirs<-temp_dirs[grep('islam|grunS|grun2i', temp_dirs)]
temp_file<-paste(temp_dirs, list.files(temp_dirs, pattern='ccc.RData'), sep='/')
names(temp_file)<-sapply(strsplit(temp_dirs, '_'), '[', 4)
for(i in 1:length(temp_file)){
  print(names(temp_file)[i])
  load(temp_file[[i]])
  nCountsMmus_list[[names(temp_file)[i]]]<-as.data.frame(t(scale(as.matrix(t(varinfo$mat)), scale=F, center=(-1)*varinfo$avmodes[rownames(varinfo$mat)])))
  nCountsMmus_mean_list[[names(temp_file)[i]]]<-varinfo$avmodes
  #normalise read counts (after cleaning counts)
  temp_sf<-estimateSizeFactorsForMatrix(cd[,rownames(knn)])
  ncd_list[[names(temp_file)[i]]]<-t( t(cd) / temp_sf )
  meansMmus_list[[names(temp_file)[i]]]<-varinfo$avmodes
  genevarrank_list[[names(temp_file)[i]]]<-genevar_rank
}

###########################
### Load and format feature data for integrated linear models
###########################

#ChIP enrichment signals
signal_type<-'_enrich_'
temp_nnfreq<-do.call("cbind", lapply(paste(chip_sub_folder, list.files(chip_sub_folder, pattern='txt'), sep="/"), read.table, header=F, sep='\t', stringsAsFactors=F))
colnames(temp_nnfreq)<-gsub(".txt", "", gsub(paste('txt', signal_type, sep=''), "", list.files(chip_sub_folder, pattern='txt')))
#Subset
temp_nnfreq<-temp_nnfreq[,-which(sapply(strsplit(colnames(temp_nnfreq), '_enrich_'), '[', 2) %in% c("GSE31039_LICR_ChipSeq_ES-Bruce4_H3K27ac_E0",
                                                                                                    "GSE31039_LICR_ChipSeq_ES-Bruce4_H3K36me3_E0",
                                                                                                    "GSE31039_LICR_ChipSeq_ES-Bruce4_H3K4me1",
                                                                                                    "GSE31039_LICR_ChipSeq_ES-Bruce4_H3K4me3",
                                                                                                    "GSE31039_LICR_ChipSeq_ES-Bruce4_H3K9ac_E0",
                                                                                                    "GSE31039_LICR_ChipSeq_ES-Bruce4_H3K9me3_E0",
                                                                                                    "GSE32218_Stanford_ChipSeq_ES-E14_H3K4me1_std",
                                                                                                    "GSE32218_Stanford_ChipSeq_ES-E14_H3K4me3_std"))]
#Simplify column names
colnames(temp_nnfreq)<-lapply(lapply(strsplit(colnames(temp_nnfreq), '_'), '[', c(1,8)), paste, collapse='_')
#Gene names
rownames(temp_nnfreq)<-names(allgeneInfo_gr)

#Maximum TATA scores for all core promoters
load(TATAscore_file)
FPR_10<-quantile(coreProm_TATAscores[names(corePromInfo_gr)[which(corePromInfo_gr$"sequence_TATAbox"==0)]], 0.9)

#All cLADs
lad<-read.table(input_file_cLAD, header=F, stringsAsFactors=F, sep='\t')
lad_gr<-GRanges(seqnames=as(lad[,1], "Rle"), ranges=IRanges(lad[,2], lad[,3]), seqlengths=seqlengths(allgeneInfo_gr))

#SE target score correlation
load(SEcor_file)

#Get Ensembl gene annotations
temp<-ensemblAPI_genes_GRanges(ensapi_file, c(1:19, 'X', 'Y', 'MT'))
allgenes<-temp[['allgenes']]
#MGI symbol to Ensembl ID mapping
genes<-unique(allgenes[,c('mgi_symbol', 'ensembl_gene_id', 'chromosome_name')])
genes<-genes[genes$chromosome_name %in% c(1:19, 'X', 'Y', 'MT'),]
#Remove ambiguous mappings (force 1:1)
genes_dup1<-genes[duplicated(genes$mgi_symbol),]$mgi_symbol
genes_dup2<-genes[duplicated(genes$ensembl_gene_id),]$ensembl_gene_id
genes<-genes[!genes$mgi_symbol %in% genes_dup1 & !genes$ensembl_gene_id %in% genes_dup2,]
#Convert to list
rownames(genes)<-genes[,1]
geneID_list<-as.list(t(genes)[2,])

#Friedel mRNA stability
decay_friedel<-read.table(mrna_decay_path_friedel, header=T, sep='\t', stringsAsFactors=F)
#Genes with Ensembl ids
decay_friedel<-decay_friedel[decay_friedel$Gene.symbol %in% names(geneID_list),]
rownames(decay_friedel)<-geneID_list[decay_friedel$Gene.symbol]
mRNA_decay<-log10(decay_friedel[,'Half.life..h.'])
names(mRNA_decay)<-rownames(decay_friedel)

#Promoter width
#Broad and sharp promoters from CAGE data (whole mouse embryo E11 - E18)
TSSclust<-read.table(cage_consensus_clusters_path, header=F, sep='\t', stringsAsFactors=F)
#Subset to canonical chromosomes
TSSclust<-TSSclust[TSSclust[,1] %in% seqnames(allgeneInfo_gr),]
#GRanges
TSSclust_gr<-trim(GRanges(seqnames=as(TSSclust[,1], "Rle"), ranges=IRanges(TSSclust[,2], TSSclust[,3]), strand=as(TSSclust[,4], "Rle"), seqlengths=seqlengths(allgeneInfo_gr)))
#Overlap with core promoters
temp_ol<-as.data.frame(findOverlaps(corePromInfo_gr, TSSclust_gr))
temp_ol$distance<-distance(tssInfo_gr[temp_ol[,1],], TSSclust_gr[temp_ol[,2],])
#Choose cluster that's closest to TSS
temp_ol<-temp_ol[order(temp_ol$distance, decreasing=F),]
temp_ol<-temp_ol[!duplicated(temp_ol[,1]),]
TSSclust_width<-rep(NA, length(corePromInfo_gr))
names(TSSclust_width)<-names(corePromInfo_gr)
TSSclust_width[temp_ol[,1]]<-width(TSSclust_gr[temp_ol[,2],])
#Mask chrM and chrY
TSSclust_width[as.character(seqnames(corePromInfo_gr)) %in% c('chrY', 'chrM')]<-NA

###########################
### Integrated models with core promoter features including promoter width
###########################

#GLMNET
set.seed(1)
feature_levels<-plot_coefficients(model_wrapper('islam', model_type='glmnet', feature_subset='coreprom_promwidth'),
                                  model_wrapper('grunS', model_type='glmnet', feature_subset='coreprom_promwidth'),
                                  model_wrapper('grun2i', model_type='glmnet', feature_subset='coreprom_promwidth'),
                                  output_file='glm_coefficients_all_coreprom_promwidth_glmnet.pdf')

#GLM_UNI
plot_coefficients(model_wrapper('islam', model_type='glm_uni', feature_subset='coreprom_promwidth'),
                  model_wrapper('grunS', model_type='glm_uni', feature_subset='coreprom_promwidth'),
                  model_wrapper('grun2i', model_type='glm_uni', feature_subset='coreprom_promwidth'),
                  output_file='glm_coefficients_all_coreprom_promwidth_glm_uni.pdf', variable_levels=feature_levels)

#GLM
plot_coefficients(model_wrapper('islam', model_type='glm', feature_subset='coreprom_promwidth'),
                  model_wrapper('grunS', model_type='glm', feature_subset='coreprom_promwidth'),
                  model_wrapper('grun2i', model_type='glm', feature_subset='coreprom_promwidth'),
                  output_file='glm_coefficients_all_coreprom_promwidth_glm.pdf', variable_levels=feature_levels)

###########################
### Integrated models with core promoter + chromatin features
###########################

#GLMNET
set.seed(1)
feature_levels<-plot_coefficients(model_wrapper('islam', model_type='glmnet', feature_subset='coreprom_chrom'),
                                  model_wrapper('grunS', model_type='glmnet', feature_subset='coreprom_chrom'),
                                  model_wrapper('grun2i', model_type='glmnet', feature_subset='coreprom_chrom'),
                                  output_file='glm_coefficients_all_coreprom_chrom_glmnet.pdf', width=10)

#GLM_UNI
plot_coefficients(model_wrapper('islam', model_type='glm_uni', feature_subset='coreprom_chrom'),
                  model_wrapper('grunS', model_type='glm_uni', feature_subset='coreprom_chrom'),
                  model_wrapper('grun2i', model_type='glm_uni', feature_subset='coreprom_chrom'),
                  output_file='glm_coefficients_all_coreprom_chrom_glm_uni.pdf', variable_levels=feature_levels, width=10)

#GLM
plot_coefficients(model_wrapper('islam', model_type='glm', feature_subset='coreprom_chrom'),
                  model_wrapper('grunS', model_type='glm', feature_subset='coreprom_chrom'),
                  model_wrapper('grun2i', model_type='glm', feature_subset='coreprom_chrom'),
                  output_file='glm_coefficients_all_coreprom_chrom_glm.pdf', variable_levels=feature_levels, width=10)

###########################
### Integrated models with core promoter + + chromatin features + chrX + mRNA stability + SE
###########################

#GLMNET
set.seed(1)
feature_levels<-plot_coefficients(model_wrapper('islam', model_type='glmnet', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                                  model_wrapper('grunS', model_type='glmnet', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                                  model_wrapper('grun2i', model_type='glmnet', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                                  output_file='glm_coefficients_all_coreprom_chrom_chrX_mRNAstability_SE_glmnet.pdf', width=10)

#GLM_UNI
plot_coefficients(model_wrapper('islam', model_type='glm_uni', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                  model_wrapper('grunS', model_type='glm_uni', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                  model_wrapper('grun2i', model_type='glm_uni', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                  output_file='glm_coefficients_all_coreprom_chrom_chrX_mRNAstability_SE_glm_uni.pdf', variable_levels=feature_levels, width=10)

#GLM
plot_coefficients(model_wrapper('islam', model_type='glm', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                  model_wrapper('grunS', model_type='glm', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                  model_wrapper('grun2i', model_type='glm', feature_subset='coreprom_chrom_X_mRNAstability_SE'),
                  output_file='glm_coefficients_all_coreprom_chrom_chrX_mRNAstability_SE_glm.pdf', variable_levels=feature_levels, width=10)
