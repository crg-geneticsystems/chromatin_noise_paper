
###########################
### REQUIRED PACKAGES
###########################

library(scde)
library(biomaRt)
library(ShortRead)
library(RColorBrewer)
library(coin)
library(ggplot2)
library(reshape2)
library(plyr)
library(GGally)
library(rtracklayer)
library(ggExtra)

###########################
### FUNCTIONS
###########################

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family="binomial"), ...)
}

###########################
### GLOBALS
###########################

#Main paths
base_path<-getwd() #Assumes current working directory is repository base folder ("chromatin_noise_paper")
misc_functions_path<-file.path(base_path, 'misc_functions.R')
misc_data_path<-file.path(base_path, 'misc_data_files')
project_name<-'04_adjvariance_bias_plots'
working_directory<-file.path(base_path, project_name)

#Miscellaneous data files
SEcor_file<-file.path(misc_data_path, 'cor_list_SE.RData')
go_file<-file.path(misc_data_path, 'GOdata_3ontologies.RData')

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

#Load GO data
load(go_file)

#Load input data (after cell cycle correction)
gst_chr_list<-list()
gst_go_list<-list()
meansMmus_list<-list()
temp_dirs<-list.dirs('..', recursive=F)
temp_dirs<-temp_dirs[grep('pagoda_adjvariance_bias', temp_dirs)]
#Retain cell cycle corrected data
temp_dirs<-temp_dirs[grep('ccc', temp_dirs)]
# temp_dirs<-temp_dirs[-grep("covariancecontrol|graf|klein|zeisel", temp_dirs)]
temp_file<-paste(temp_dirs, list.files(temp_dirs, pattern='ccc.RData'), sep='/')
names(temp_file)<-sapply(strsplit(temp_dirs, '_'), '[', 4)
for(i in 1:length(temp_file)){
  print(names(temp_file)[i])
  load(temp_file[[i]])
  ##GO
  gs_tests_all_GO[['biological_process']]<-gs_tests_all_GO[['biological_process']] 
  gs_tests_all_GO[['cellular_component']]<-gs_tests_all_GO[['cellular_component']]  
  gs_tests_all_GO[['molecular_function']]<-gs_tests_all_GO[['molecular_function']] 
  #Save
  gst_go_list[[names(temp_file)[i]]]<-gs_tests_all_GO
  ##Chromatin
  gs_tests_all<-gs_tests_all
  #Save
  gst_chr_list[[names(temp_file)[i]]]<-gs_tests_all
  meansMmus_list[[names(temp_file)[i]]]<-varinfo$avmodes
}

###########################
### Plot gene sets adjusted variance bias bubble plots - all
###########################

#All features total variability
plot_df<-do.call('rbind', gst_chr_list)
plot_df$dataset<-as.factor(rep(names(gst_chr_list), times=as.numeric(sapply(gst_chr_list, dim)[1,])))
#Exclude datasets
# plot_df<-plot_df[!plot_df$dataset %in% c('grafS1', 'kleinS1'),]
#Feature category
plot_df$feat_category<-extract_feature_category(rownames(plot_df))
#Feature name
plot_df$feat_name<-extract_feature_name(rownames(plot_df))
#Bubble colour
plot_df$location<-'1 TSS'
plot_df$location[grep('corePromInfo', rownames(plot_df))]<-'2 core promoter'
plot_df$location[grep('promInfo', rownames(plot_df))]<-'3 promoter'
plot_df$location[grep('allgeneInfo', rownames(plot_df))]<-'4 whole gene body'
plot_df$location<-as.character(plot_df$location)
plot_cols<-c(brewer.pal(5, 'Spectral')[-3])
names(plot_cols)<-c('1 TSS', '2 core promoter', '3 promoter', '4 whole gene body')
#Set non-signficant gene sets to smaller dot
plot_df$significance<-'FDR<10%'
plot_df$significance[plot_df$FDR>=0.1]<-'FDR>=10%'
plot_df$significance<-as.factor(plot_df$significance)
plot_shapes<-c(1, 19)
names(plot_shapes)<-c('FDR>=10%','FDR<10%')
#Plot index
temp<-plot_df[!duplicated(plot_df$feat_name),]
temp_order<-order(temp$feat_category, temp$feat_name)
plot_df$index<-factor(plot_df$feat_name, levels=temp$feat_name[temp_order])
#Plot
pdf('gene_set_variability_bubble_chrom_ccc.pdf', width=20, height=80)
d <- ggplot(plot_df, aes(index, effect_size)) + 
  geom_point(aes(color = location, shape=significance), size=4) +
  scale_colour_manual(values = plot_cols) + 
  scale_shape_manual(values = plot_shapes) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, linetype=2, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d + facet_wrap(~dataset, scales='free_x', nrow=8, ncol=1)
dev.off()

###########################
### Plot gene sets adjusted variance bias bubble plots - feature subset
###########################

#Feature subset total variability (no repeats)
gst_chr_list_subset<-lapply(gst_chr_list, feature_subset)
plot_df<-do.call('rbind', gst_chr_list_subset)
plot_df$dataset<-rep(names(gst_chr_list_subset), times=as.numeric(sapply(gst_chr_list_subset, dim)[1,]))
#Subset datasets
plot_df<-plot_df[plot_df$dataset %in% c('islam', 'grun2i', 'grunS'),]
plot_df$dataset[plot_df$dataset=='islam']<-'1 islam'
plot_df$dataset[plot_df$dataset=='grunS']<-'2 grunS'
plot_df$dataset[plot_df$dataset=='grun2i']<-'3 grun2i'
plot_df$dataset<-as.factor(plot_df$dataset)
#Feature category
plot_df$feat_category<-extract_feature_category(rownames(plot_df))
#Feature name
plot_df$feat_name<-extract_feature_name(rownames(plot_df))
#Bubble colour
plot_df$location<-'1 TSS'
plot_df$location[grep('corePromInfo', rownames(plot_df))]<-'2 core promoter'
plot_df$location[grep('promInfo', rownames(plot_df))]<-'3 promoter'
plot_df$location[grep('allgeneInfo', rownames(plot_df))]<-'4 whole gene body'
plot_df$location<-as.character(plot_df$location)
plot_cols<-c(brewer.pal(5, 'Spectral')[-3])
names(plot_cols)<-c('1 TSS', '2 core promoter', '3 promoter', '4 whole gene body')
#Set non-signficant gene sets to smaller dot
plot_df$significance<-'FDR<10%'
plot_df$significance[plot_df$FDR>=0.1]<-'FDR>=10%'
plot_df$significance<-as.factor(plot_df$significance)
plot_shapes<-c(1, 19)
names(plot_shapes)<-c('FDR>=10%','FDR<10%')
#Index
meds<-ddply(plot_df, c("feat_category", "feat_name", "location"), summarise, effect_size_median = median(effect_size))
meds<-meds[order(meds$feat_category, abs(meds$effect_size_median-0.5), decreasing=T),]
meds<-meds[!duplicated(meds$feat_name),]
meds<-meds[order(meds$feat_category, meds$effect_size_median, decreasing=T),]
plot_df$index<-factor(plot_df$feat_name, levels=meds$feat_name)
#Plot (UMI)
pdf('gene_set_variability_bubble_chromsubset_ccc.pdf', width=15, height=7)
d <- ggplot(plot_df, aes(index, effect_size)) + 
  geom_point(aes(color = location, shape=significance), size=4) +
  scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  scale_shape_manual(values = plot_shapes) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d + facet_wrap(~dataset, ncol=3, nrow=1)
dev.off()

###########################
### Plot gene sets adjusted variance bias dot plots - feature subset - datasetsymbols
###########################

#Feature subset total variability (no repeats)
gst_chr_list_subset<-lapply(gst_chr_list, feature_subset)
plot_df<-do.call('rbind', gst_chr_list_subset)
plot_df$dataset<-rep(names(gst_chr_list_subset), times=as.numeric(sapply(gst_chr_list_subset, dim)[1,]))
#Subset datasets
plot_df<-plot_df[plot_df$dataset %in% c('islam', 'grunS', 'grun2i'),]
#Feature category
plot_df$feat_category<-extract_feature_category(rownames(plot_df))
#Feature name
plot_df$feat_name<-extract_feature_name(rownames(plot_df))
#Bubble colour
plot_df$location<-'1 TSS'
plot_df$location[grep('corePromInfo', rownames(plot_df))]<-'2 core promoter'
plot_df$location[grep('promInfo', rownames(plot_df))]<-'3 promoter'
plot_df$location[grep('allgeneInfo', rownames(plot_df))]<-'4 whole gene body'
plot_df$location<-as.character(plot_df$location)
plot_cols<-c(brewer.pal(5, 'Spectral')[-3])
names(plot_cols)<-c('1 TSS', '2 core promoter', '3 promoter', '4 whole gene body')
#Set non-signficant gene sets to smaller dot
plot_df$significance<-'FDR<10%'
plot_df$significance[plot_df$FDR>=0.1]<-'FDR>=10%'
plot_df$significance<-as.factor(plot_df$significance)
plot_shapes<-c(1, 19)
names(plot_shapes)<-c('FDR>=10%','FDR<10%')
#Index
meds<-ddply(plot_df, c("feat_category", "feat_name", "location"), summarise, effect_size_median = median(effect_size))
meds<-meds[order(meds$feat_category, abs(meds$effect_size_median-0.5), decreasing=T),]
meds<-meds[!duplicated(meds$feat_name),]
meds<-meds[order(meds$feat_category, meds$effect_size_median, decreasing=T),]
plot_df$index<-factor(plot_df$feat_name, levels=meds$feat_name)
#Plot
pdf('gene_set_variability_points_decsubset_ccc_datasetsymbols.pdf', width=5, height=5)
d <- ggplot(plot_df[plot_df$feat_category=='Sequence',], aes(index, effect_size, group = interaction(index, location))) + 
  #   geom_boxplot(aes(color = location)) + 
  geom_point(aes(color = location, shape = dataset), size=2, position = position_dodge(width = 0.75)) +
  geom_line(aes(color = location), position = position_dodge(width = 0.75)) +
  # scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d
d <- ggplot(plot_df[plot_df$feat_category=='Local chromatin (peaks)',], aes(index, effect_size, group = interaction(index, location))) + 
  #   geom_boxplot(aes(color = location)) + 
  geom_point(aes(color = location, shape = dataset), size=2, position = position_dodge(width = 0.75)) +
  geom_line(aes(color = location), position = position_dodge(width = 0.75)) +
  # scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d
d <- ggplot(plot_df[plot_df$feat_category=='Global chromatin',], aes(index, effect_size, group = interaction(index, location))) + 
  #   geom_boxplot(aes(color = location)) + 
  geom_point(aes(color = location, shape = dataset), size=2, position = position_dodge(width = 0.75)) +
  geom_line(aes(color = location), position = position_dodge(width = 0.75)) +
  # scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d
dev.off()
#Save order
meds_noise<-meds
#Save plot_df
plot_df_real<-plot_df

###########################
### Plot gene sets adjusted variance bias dot plots - feature subset - datasetsymbols (pool-and-split)
###########################

#Feature subset total variability (no repeats)
gst_chr_list_subset<-lapply(gst_chr_list, feature_subset)
plot_df<-do.call('rbind', gst_chr_list_subset)
plot_df$dataset<-rep(names(gst_chr_list_subset), times=as.numeric(sapply(gst_chr_list_subset, dim)[1,]))
#Subset datasets
plot_df<-plot_df[plot_df$dataset %in% c('grunSPandS', 'grun2iPandS'),]
#Feature category
plot_df$feat_category<-extract_feature_category(rownames(plot_df))
#Feature name
plot_df$feat_name<-extract_feature_name(rownames(plot_df))
#Bubble colour
plot_df$location<-'1 TSS'
plot_df$location[grep('corePromInfo', rownames(plot_df))]<-'2 core promoter'
plot_df$location[grep('promInfo', rownames(plot_df))]<-'3 promoter'
plot_df$location[grep('allgeneInfo', rownames(plot_df))]<-'4 whole gene body'
plot_df$location<-as.character(plot_df$location)
plot_cols<-c(brewer.pal(5, 'Spectral')[-3])
names(plot_cols)<-c('1 TSS', '2 core promoter', '3 promoter', '4 whole gene body')
#Set non-signficant gene sets to smaller dot
plot_df$significance<-'FDR<10%'
plot_df$significance[plot_df$FDR>=0.1]<-'FDR>=10%'
plot_df$significance<-as.factor(plot_df$significance)
plot_shapes<-c(1, 19)
names(plot_shapes)<-c('FDR>=10%','FDR<10%')
#Index
meds<-meds_noise
plot_df$index<-factor(plot_df$feat_name, levels=meds$feat_name)
#Plot
pdf('gene_set_variability_points_decsubset_ccc_fixedscale2_datasetsymbols_PandS.pdf', width=5, height=5)
d <- ggplot(plot_df[plot_df$feat_category=='Sequence',], aes(index, effect_size, group = interaction(index, location))) + 
  #   geom_boxplot(aes(color = location)) + 
  geom_point(data = plot_df_real[plot_df$feat_category=='Sequence',], colour=NA) +
  geom_point(aes(color = location, shape = dataset), size=2, position = position_dodge(width = 0.75)) +
  geom_line(aes(color = location), position = position_dodge(width = 0.75)) +
  # scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d
d <- ggplot(plot_df[plot_df$feat_category=='Local chromatin (peaks)',], aes(index, effect_size, group = interaction(index, location))) + 
  #   geom_boxplot(aes(color = location)) + 
  geom_point(data = plot_df_real[plot_df$feat_category=='Local chromatin (peaks)',], colour=NA) +
  geom_point(aes(color = location, shape = dataset), size=2, position = position_dodge(width = 0.75)) +
  geom_line(aes(color = location), position = position_dodge(width = 0.75)) +
  # scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d
d <- ggplot(plot_df[plot_df$feat_category=='Global chromatin',], aes(index, effect_size, group = interaction(index, location))) + 
  #   geom_boxplot(aes(color = location)) + 
  geom_point(data = plot_df_real[plot_df$feat_category=='Global chromatin',], colour=NA) +
  geom_point(aes(color = location, shape = dataset), size=2, position = position_dodge(width = 0.75)) +
  geom_line(aes(color = location), position = position_dodge(width = 0.75)) +
  # scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d
dev.off()

###########################
### Plot gene sets adjusted variance bias bubble plots - feature subset (pool-and-split)
###########################

#Feature subset total variability (no repeats)
gst_chr_list_subset<-lapply(gst_chr_list, feature_subset)
plot_df<-do.call('rbind', gst_chr_list_subset)
plot_df$dataset<-rep(names(gst_chr_list_subset), times=as.numeric(sapply(gst_chr_list_subset, dim)[1,]))
#Subset datasets
plot_df<-plot_df[plot_df$dataset %in% c('islam', 'grun2iPandS', 'grunSPandS'),]
plot_df$dataset[plot_df$dataset=='islam']<-'1 islam'
plot_df$dataset[plot_df$dataset=='grunSPandS']<-'2 grunSPandS'
plot_df$dataset[plot_df$dataset=='grun2iPandS']<-'3 grun2iPandS'
plot_df$dataset<-as.factor(plot_df$dataset)
#Feature category
plot_df$feat_category<-extract_feature_category(rownames(plot_df))
#Feature name
plot_df$feat_name<-extract_feature_name(rownames(plot_df))
#Bubble colour
plot_df$location<-'1 TSS'
plot_df$location[grep('corePromInfo', rownames(plot_df))]<-'2 core promoter'
plot_df$location[grep('promInfo', rownames(plot_df))]<-'3 promoter'
plot_df$location[grep('allgeneInfo', rownames(plot_df))]<-'4 whole gene body'
plot_df$location<-as.character(plot_df$location)
plot_cols<-c(brewer.pal(5, 'Spectral')[-3])
names(plot_cols)<-c('1 TSS', '2 core promoter', '3 promoter', '4 whole gene body')
#Set non-signficant gene sets to smaller dot
plot_df$significance<-'FDR<10%'
plot_df$significance[plot_df$FDR>=0.1]<-'FDR>=10%'
plot_df$significance<-as.factor(plot_df$significance)
plot_shapes<-c(1, 19)
names(plot_shapes)<-c('FDR>=10%','FDR<10%')
#Plot index
plot_df$index<-factor(plot_df$feat_name, levels=meds$feat_name)
#Plot (UMI)
pdf('gene_set_variability_bubble_chromsubset_ccc_PandS.pdf', width=15, height=7)
d <- ggplot(plot_df, aes(index, effect_size)) + 
  geom_point(aes(color = location, shape=significance), size=4) +
  scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  scale_shape_manual(values = plot_shapes) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d + facet_wrap(~dataset, ncol=3, nrow=1)
dev.off()

###########################
### Plot smoothed adjusted variance for selected chromatin gene sets
###########################

#Create output directory
dir.create("smoothed_adjvariance_plots", showWarnings = FALSE)

#Adjusted variance plots of selected trends
#Load adjusted variance data for selected gene sets (after cell cycle correction)
plot_list<-list()
temp_dirs<-list.dirs('..', recursive=F)
temp_dirs<-temp_dirs[grep('pagoda_adjvariance_bias', temp_dirs)]
#Retain cell cycle corrected data
temp_dirs<-temp_dirs[grep('ccc', temp_dirs)]
# temp_dirs<-temp_dirs[-grep("covariancecontrol|graf|klein|zeisel", temp_dirs)]
#Datasets of interest
temp_dirs<-temp_dirs[grep('islam|grun', temp_dirs)]
temp_file<-paste(temp_dirs, list.files(temp_dirs, pattern='ccc.RData'), sep='/')
names(temp_file)<-sapply(strsplit(temp_dirs, '_'), '[', 4)
for(i in 1:length(temp_file)){
  print(names(temp_file)[i])
  load(temp_file[[i]])
  plot_list[[names(temp_file)[i]]]<-data.frame(mean_expr=varinfo$avmodes, 
                                               adjusted_variance=genevar_rank[names(varinfo$avmodes)], 
                                               TATA_coreprom=as.factor(corePromInfo_gr[names(varinfo$avmodes),]$"sequence_TATAbox"),
                                               CGI_coreprom=as.factor(corePromInfo_gr[names(varinfo$avmodes),]$"sequence_CpG island"),
                                               SharpTSS_coreprom=as.factor(corePromInfo_gr[names(varinfo$avmodes),]$"sequence_Sharp TSS"),
                                               BroadTSS_coreprom=as.factor(corePromInfo_gr[names(varinfo$avmodes),]$"sequence_Broad TSS"),
                                               H3K27me3_tss=as.factor(tssInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k27me3"),
                                               H3K27me3only_prom=as.factor(as.numeric(promInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k27me3" & !promInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k04me3")),
                                               Bivalent_prom=as.factor(as.numeric(promInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k27me3" & promInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k04me3")),
                                               H3K9me3_prom=as.factor(promInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k09me3"),
                                               H3K4me3only_prom=as.factor(as.numeric(!promInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k27me3" & promInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k04me3")),
                                               H3K36me3_genebody=as.factor(allgeneInfo_gr[names(varinfo$avmodes),]$"localchromatin_peak_H3k36me3"))
}

#All UMI datasets
plot_df<-do.call('rbind', plot_list)
plot_df$dataset<-rep(names(plot_list), times=as.numeric(sapply(plot_list, dim)[1,]))
#Subset datasets
plot_df<-plot_df[plot_df$dataset %in% c('islam', 'grun2i', 'grunS'),]
plot_df$dataset[plot_df$dataset=='islam']<-'1 islam'
plot_df$dataset[plot_df$dataset=='grunS']<-'2 grunS'
plot_df$dataset[plot_df$dataset=='grun2i']<-'3 grun2i'
plot_df$dataset<-as.factor(plot_df$dataset)
plot_df$adjusted_variance<-(plot_df$adjusted_variance-1)/100
#Plot
for(i in colnames(plot_list[[1]])[3:dim(plot_list[[1]])[2]]){
  temp_col<-c(2,2,2,2,1,4,4,4,4,5)
  names(temp_col)<-c('TATA_coreprom', 'CGI_coreprom', 'SharpTSS_coreprom', 'BroadTSS_coreprom', 'H3K27me3_tss', 'H3K27me3only_prom', 'Bivalent_prom', 'H3K9me3_prom', 'H3K4me3only_prom', 'H3K36me3_genebody')
  plot_cols<-c(brewer.pal(5, 'Spectral')[temp_col[i]], 'black')
  names(plot_cols)<-c('1', '0')
  print(i)
  temp_plot_df<-plot_df[,c('mean_expr', 'adjusted_variance', i, 'dataset')]
  colnames(temp_plot_df)[3]<-'gene_set'
  d <- ggplot(temp_plot_df[!is.na(temp_plot_df$gene_set),], aes(mean_expr, adjusted_variance, colour=gene_set)) + scale_x_log10() + 
    binomial_smooth(formula = y ~ splines::ns(x, 2), se=T, size=1, aes(fill = gene_set)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = plot_cols) + 
    scale_fill_manual(values = plot_cols) + 
    geom_rug(data=temp_plot_df[which(temp_plot_df$gene_set==1),], col=adjustcolor(plot_cols[1],alpha.f=0.3), sides='b') +
    theme_bw()
  d + facet_wrap(~dataset, ncol=3, nrow=1)
  ggsave(file=paste('smoothed_adjvariance_plots/gene_set_adjvariance_UMI_', i, '.pdf', sep=''), width=15, height=5)
}

#All UMI datasets (pool-and-split controls)
plot_df<-do.call('rbind', plot_list)
plot_df$dataset<-rep(names(plot_list), times=as.numeric(sapply(plot_list, dim)[1,]))
#Subset datasets
plot_df<-plot_df[plot_df$dataset %in% c('islam', 'grun2iPandS', 'grunSPandS'),]
plot_df$dataset[plot_df$dataset=='islam']<-'1 islam'
plot_df$dataset[plot_df$dataset=='grunSPandS']<-'2 grunSPandS'
plot_df$dataset[plot_df$dataset=='grun2iPandS']<-'3 grun2iPandS'
plot_df$dataset<-as.factor(plot_df$dataset)
plot_df$adjusted_variance<-(plot_df$adjusted_variance-1)/100
#Plot
for(i in colnames(plot_list[[1]])[3:10]){
  temp_col<-c(2,2,2,2,1,4,4,4,4,5)
  names(temp_col)<-c('TATA_coreprom', 'CGI_coreprom', 'SharpTSS_coreprom', 'BroadTSS_coreprom', 'H3K27me3_tss', 'H3K27me3only_prom', 'Bivalent_prom', 'H3K9me3_prom', 'H3K4me3only_prom', 'H3K36me3_genebody')
  plot_cols<-c(brewer.pal(5, 'Spectral')[temp_col[i]], 'black')
  names(plot_cols)<-c('1', '0')
  print(i)
  temp_plot_df<-plot_df[,c('mean_expr', 'adjusted_variance', i, 'dataset')]
  colnames(temp_plot_df)[3]<-'gene_set'
  d <- ggplot(temp_plot_df[!is.na(temp_plot_df$gene_set),], aes(mean_expr, adjusted_variance, colour=gene_set)) + scale_x_log10() + 
    binomial_smooth(formula = y ~ splines::ns(x, 2), se=T, size=1, aes(fill = gene_set)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = plot_cols) + 
    scale_fill_manual(values = plot_cols) + 
    theme_bw()
  d + facet_wrap(~dataset, ncol=3, nrow=1)
  ggsave(file=paste('smoothed_adjvariance_plots/gene_set_adjvariance_UMIPandS_', i, '.pdf', sep=''), width=15, height=5)
}

###########################
### Plot smoothed adjusted variance for mismatched chromatin gene sets
###########################

load('../pagoda_adjvariance_bias_islam_ccc/pagoda_adjvariance_bias_islam_ccc.RData')
mismatch_histone_df<-data.frame(H3K27me3=promInfo_gr$"localchromatin_peak_H3k27me3"==1, #repressive
                                H3K4me1=tssInfo_gr$"localchromatin_peak_H3k04me1"==1, #repressive
                                H3K9me3=tssInfo_gr$"localchromatin_peak_H3k09me3"==1, #repressive
                                H3K27ac=promInfo_gr$"localchromatin_peak_H3k27ac"==0, #active
                                H3K9ac=allgeneInfo_gr$"localchromatin_peak_H3k09ac"==0, #active
                                H3K36me3=allgeneInfo_gr$"localchromatin_peak_H3k36me3"==0, #active
                                H3K4me3=allgeneInfo_gr$"localchromatin_peak_H3k04me3"==0) #active
rownames(mismatch_histone_df)<-names(allgeneInfo_gr)
mismatch_histone_count<-apply(mismatch_histone_df, 1, sum, na.rm=T)

#Adjusted variance plots of quantitatively mismatched chromatin (combinatorial)
#Load adjusted variance data for selected gene sets (after cell cycle correction)
plot_list<-list()
temp_dirs<-list.dirs('..', recursive=F)
temp_dirs<-temp_dirs[grep('pagoda_adjvariance_bias', temp_dirs)]
#Retain cell cycle corrected data
temp_dirs<-temp_dirs[grep('ccc', temp_dirs)]
# temp_dirs<-temp_dirs[-grep("covariancecontrol|graf|klein|zeisel", temp_dirs)]
#Datasets of interest
temp_dirs<-temp_dirs[grep('islam|grun', temp_dirs)]
temp_file<-paste(temp_dirs, list.files(temp_dirs, pattern='ccc.RData'), sep='/')
names(temp_file)<-sapply(strsplit(temp_dirs, '_'), '[', 4)
for(i in 1:length(temp_file)){
  print(names(temp_file)[i])
  load(temp_file[[i]])
  plot_list[[names(temp_file)[i]]]<-data.frame(mean_expr=varinfo$avmodes, 
                                               adjusted_variance=genevar_rank[names(varinfo$avmodes)],
                                               mismatch_histone_count=mismatch_histone_count[names(varinfo$avmodes)],
                                               cgi_promoter=corePromInfo_gr[names(varinfo$avmodes),]$"sequence_CpG island")
}

#All datasets
plot_df<-do.call('rbind', plot_list)
plot_df$dataset<-rep(names(plot_list), times=as.numeric(sapply(plot_list, dim)[1,]))
plot_df<-plot_df[plot_df$dataset %in% c('islam', 'grun2i', 'grunS'),]
plot_df$dataset[plot_df$dataset=='islam']<-'1 islam'
plot_df$dataset[plot_df$dataset=='grunS']<-'2 grunS'
plot_df$dataset[plot_df$dataset=='grun2i']<-'3 grun2i'
plot_df$dataset<-as.factor(plot_df$dataset)
plot_df$mismatch_histone_count[plot_df$mismatch_histone_count>=3]<-'3+'
plot_df$mismatch_histone_count<-as.factor(plot_df$mismatch_histone_count)
plot_df$adjusted_variance<-(plot_df$adjusted_variance-1)/100
plot_cols<-c(rev(brewer.pal(5, 'Spectral')[-3]))
names(plot_cols)<-c(0:2, '3+')
#Plot - histone
d <- ggplot(plot_df, aes(mean_expr, adjusted_variance, colour=mismatch_histone_count)) + scale_x_log10() + 
  binomial_smooth(formula = y ~ splines::ns(x, 2), se=T, size=1, aes(fill = mismatch_histone_count), alpha = 0.2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = plot_cols) + 
  scale_fill_manual(values = plot_cols) + 
  theme_bw()
d + facet_wrap(~dataset, ncol=3, nrow=1)
ggsave(file=paste('smoothed_adjvariance_plots/gene_set_adjvariance_histone_mismatchcount.pdf', sep=''), width=15, height=5)
#Plot - CGI promoters
d <- ggplot(plot_df[which(plot_df$cgi_promoter==1),], aes(mean_expr, adjusted_variance, colour=mismatch_histone_count)) + scale_x_log10() + 
  binomial_smooth(formula = y ~ splines::ns(x, 2), se=T, size=1, aes(fill = mismatch_histone_count), alpha = 0.2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = plot_cols) + 
  scale_fill_manual(values = plot_cols) + 
  theme_bw()
d + facet_wrap(~dataset, ncol=3, nrow=1)
ggsave(file=paste('smoothed_adjvariance_plots/gene_set_adjvariance_histone_mismatchcount_cgiprom.pdf', sep=''), width=15, height=5)
#Plot - non CGI promoters
d <- ggplot(plot_df[which(plot_df$cgi_promoter==0),], aes(mean_expr, adjusted_variance, colour=mismatch_histone_count)) + scale_x_log10() + 
  binomial_smooth(formula = y ~ splines::ns(x, 2), se=T, size=1, aes(fill = mismatch_histone_count), alpha = 0.2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = plot_cols) + 
  scale_fill_manual(values = plot_cols) + 
  theme_bw()
d + facet_wrap(~dataset, ncol=3, nrow=1)
ggsave(file=paste('smoothed_adjvariance_plots/gene_set_adjvariance_histone_mismatchcount_noncgiprom.pdf', sep=''), width=15, height=5)

