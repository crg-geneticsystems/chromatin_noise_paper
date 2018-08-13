
###########################
### REQUIRED PACKAGES
###########################

library(optparse)
library(scde)
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
misc_data_path<-file.path(base_path, 'misc_data_files')
project_name<-'06_mRNAstability_plots'
working_directory<-file.path(base_path, project_name)

#Gene annotations and gene set files
ensapi_file<-file.path(misc_data_path, 'ens71_allgenefeatures.txt')
ensapi_refseq_file<-file.path(misc_data_path, 'ens71_allgenefeatures_refseq.txt')
geneset_file<-file.path(misc_data_path, 'SC_preprocessing_genesets.RData')

#mRNA decay files
mrna_decay_path_sharova<-file.path(misc_data_path, 'DNA_Res_2009_Sharova.txt')
mrna_decay_path_friedel<-file.path(misc_data_path, 'Friedel2009_NIH3T3_mRNAhalflife_unix.txt')
mrna_decay_path_tippmann<-file.path(misc_data_path, 'ES_TN_halflives.txt')
mrna_decay_path_schwanhausser<-file.path(misc_data_path, 'Nature_2011_Schwanhausser_unix.txt')

#Colour palette
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
#Set working directory
setwd(working_directory)

#Required functions
source(misc_functions_path)

###########################
### Load required data and annotations
###########################

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

#Load gene set data
load(geneset_file)
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
### Test chromatin gene set mRNA stability bias
###########################

#Initialise list of results
gst_dec_all_list<-list()

### Sharova et al.
###########################

decay_sharova<-read.table(mrna_decay_path_sharova, header=T, sep='\t', stringsAsFactors=F, skip=5, comment.char="", quote=NULL)
#Genes with Ensembl ids
decay_sharova<-decay_sharova[decay_sharova$geneSymbol %in% names(geneID_list),]
rownames(decay_sharova)<-geneID_list[decay_sharova$geneSymbol]
#Test
gst_dec_all<-list()
gst_dec_all[['tssInfo']]<-genesets_biased(as.list(elementMetadata(tssInfo_gr[rownames(decay_sharova),])), decay_sharova$bined.1)
gst_dec_all[['corePromInfo']]<-genesets_biased(as.list(elementMetadata(corePromInfo_gr[rownames(decay_sharova),])), decay_sharova$bined.1)
gst_dec_all[['promInfo']]<-genesets_biased(as.list(elementMetadata(promInfo_gr[rownames(decay_sharova),])), decay_sharova$bined.1)
gst_dec_all[['allgeneInfo']]<-genesets_biased(as.list(elementMetadata(allgeneInfo_gr[rownames(decay_sharova),])), decay_sharova$bined.1)
gst_dec_all_list[['sharova']]<-do.call('rbind', gst_dec_all)
#Upper quartile stable mRNAs
stable_RNA_sharova<-rownames(decay_sharova)[decay_sharova$bined.1>quantile(decay_sharova$bined.1)[4]]
#Lower quartile unstable mRNAs
unstable_RNA_sharova<-rownames(decay_sharova)[decay_sharova$bined.1<quantile(decay_sharova$bined.1)[2]]

### Friedel et al.
###########################

decay_friedel<-read.table(mrna_decay_path_friedel, header=T, sep='\t', stringsAsFactors=F)
#Genes with Ensembl ids
decay_friedel<-decay_friedel[decay_friedel$Gene.symbol %in% names(geneID_list),]
rownames(decay_friedel)<-geneID_list[decay_friedel$Gene.symbol]
#Test
gst_dec_all<-list()
gst_dec_all[['tssInfo']]<-genesets_biased(as.list(elementMetadata(tssInfo_gr[rownames(decay_friedel),])), decay_friedel$Half.life..h.)
gst_dec_all[['corePromInfo']]<-genesets_biased(as.list(elementMetadata(corePromInfo_gr[rownames(decay_friedel),])), decay_friedel$Half.life..h.)
gst_dec_all[['promInfo']]<-genesets_biased(as.list(elementMetadata(promInfo_gr[rownames(decay_friedel),])), decay_friedel$Half.life..h.)
gst_dec_all[['allgeneInfo']]<-genesets_biased(as.list(elementMetadata(allgeneInfo_gr[rownames(decay_friedel),])), decay_friedel$Half.life..h.)
gst_dec_all_list[['friedel']]<-do.call('rbind', gst_dec_all)

### Tippmann et al.
###########################

decay_tippmann<-read.table(mrna_decay_path_tippmann, header=T, sep=' ', stringsAsFactors=F)
ens71_refseq<-read.table(ensapi_refseq_file, header=F, sep='\t', stringsAsFactors=F)
ens71_refseq[,8]<-sapply(strsplit(ens71_refseq[,8], '\\.'), '[', 1)
ens71_refseq<-unique(ens71_refseq[ens71_refseq[,8]!='' & ens71_refseq[,1] %in% names(allgeneInfo_gr) & ens71_refseq[,8] %in% decay_tippmann[,1],c(1,8)])
rownames(ens71_refseq)<-ens71_refseq[,2]
#Genes with Ensembl ids
decay_tippmann<-decay_tippmann[decay_tippmann[,1] %in% rownames(ens71_refseq),]
rownames(decay_tippmann)<-ens71_refseq[decay_tippmann[,1],1]
#Test
gst_dec_all<-list()
gst_dec_all[['tssInfo']]<-genesets_biased(as.list(elementMetadata(tssInfo_gr[rownames(decay_tippmann),])), decay_tippmann$HL_ES)
gst_dec_all[['corePromInfo']]<-genesets_biased(as.list(elementMetadata(corePromInfo_gr[rownames(decay_tippmann),])), decay_tippmann$HL_ES)
gst_dec_all[['promInfo']]<-genesets_biased(as.list(elementMetadata(promInfo_gr[rownames(decay_tippmann),])), decay_tippmann$HL_ES)
gst_dec_all[['allgeneInfo']]<-genesets_biased(as.list(elementMetadata(allgeneInfo_gr[rownames(decay_tippmann),])), decay_tippmann$HL_ES)
gst_dec_all_list[['tippmann']]<-do.call('rbind', gst_dec_all)
#Upper quartile stable mRNAs
stable_RNA_tippmann<-rownames(decay_tippmann)[decay_tippmann$HL_ES>quantile(decay_tippmann$HL_ES)[4]]
#Lower quartile unstable mRNAs
unstable_RNA_tippmann<-rownames(decay_tippmann)[decay_tippmann$HL_ES<quantile(decay_tippmann$HL_ES)[2]]

### Schwanhausser et al.
###########################

decay_schwanhausser<-read.table(mrna_decay_path_schwanhausser, header=T, sep='\t', stringsAsFactors=F, comment.char="", quote=NULL)
ens71_refseq<-read.table(ensapi_refseq_file, header=F, sep='\t', stringsAsFactors=F)
ens71_refseq[,8]<-sapply(strsplit(ens71_refseq[,8], '\\.'), '[', 1)
ens71_refseq<-unique(ens71_refseq[ens71_refseq[,8]!='' & ens71_refseq[,1] %in% names(allgeneInfo_gr) & ens71_refseq[,8] %in% decay_schwanhausser$Refseq.mRNA.ID,c(1,8)])
ens71_refseq<-ens71_refseq[!ens71_refseq[,1] %in% ens71_refseq[,1][duplicated(ens71_refseq[,1])],]
rownames(ens71_refseq)<-ens71_refseq[,2]
#Genes with Ensembl ids
decay_schwanhausser<-decay_schwanhausser[decay_schwanhausser$Refseq.mRNA.ID %in% rownames(ens71_refseq),]
rownames(decay_schwanhausser)<-ens71_refseq[decay_schwanhausser$Refseq.mRNA.ID,1]
#Test
gst_dec_all<-list()
gst_dec_all[['tssInfo']]<-genesets_biased(as.list(elementMetadata(tssInfo_gr[rownames(decay_schwanhausser),])), decay_schwanhausser$mRNA.half.life.average..h.)
gst_dec_all[['corePromInfo']]<-genesets_biased(as.list(elementMetadata(corePromInfo_gr[rownames(decay_schwanhausser),])), decay_schwanhausser$mRNA.half.life.average..h.)
gst_dec_all[['promInfo']]<-genesets_biased(as.list(elementMetadata(promInfo_gr[rownames(decay_schwanhausser),])), decay_schwanhausser$mRNA.half.life.average..h.)
gst_dec_all[['allgeneInfo']]<-genesets_biased(as.list(elementMetadata(allgeneInfo_gr[rownames(decay_schwanhausser),])), decay_schwanhausser$mRNA.half.life.average..h.)
gst_dec_all_list[['schwanhausser']]<-do.call('rbind', gst_dec_all)
#Upper quartile stable mRNAs
stable_RNA_schwanhausser<-rownames(decay_schwanhausser)[decay_schwanhausser$mRNA.half.life.average..h.>quantile(decay_schwanhausser$mRNA.half.life.average..h.)[4]]
#Lower quartile unstable mRNAs
unstable_RNA_schwanhausser<-rownames(decay_schwanhausser)[decay_schwanhausser$mRNA.half.life.average..h.<quantile(decay_schwanhausser$mRNA.half.life.average..h.)[2]]

#Friedel dataset shows strongest association between stability and occurrence on the X-chromosome
temp<-do.call('rbind', gst_dec_all_list)
temp[grep('X-chromosome', rownames(temp)),]

###########################
### Test mRNA stability gene set noise bias
###########################

#Initialise object
allgeneInfo_decay_gr<-allgeneInfo_gr
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA stableSh']<-as.numeric(names(allgeneInfo_decay_gr) %in% stable_RNA_sharova)
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA stableSh_background']<-as.numeric(names(allgeneInfo_decay_gr) %in% rownames(decay_sharova))
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA unstableSh']<-as.numeric(names(allgeneInfo_decay_gr) %in% unstable_RNA_sharova)
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA unstableSh_background']<-as.numeric(names(allgeneInfo_decay_gr) %in% rownames(decay_sharova))
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA stableT']<-as.numeric(names(allgeneInfo_decay_gr) %in% stable_RNA_tippmann)
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA stableT_background']<-as.numeric(names(allgeneInfo_decay_gr) %in% rownames(decay_tippmann))
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA unstableT']<-as.numeric(names(allgeneInfo_decay_gr) %in% unstable_RNA_tippmann)
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA unstableT_background']<-as.numeric(names(allgeneInfo_decay_gr) %in% rownames(decay_tippmann))
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA stableS']<-as.numeric(names(allgeneInfo_decay_gr) %in% stable_RNA_schwanhausser)
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA stableS_background']<-as.numeric(names(allgeneInfo_decay_gr) %in% rownames(decay_schwanhausser))
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA unstableS']<-as.numeric(names(allgeneInfo_decay_gr) %in% unstable_RNA_schwanhausser)
elementMetadata(allgeneInfo_decay_gr)[,'sequence_mRNA unstableS_background']<-as.numeric(names(allgeneInfo_decay_gr) %in% rownames(decay_schwanhausser))
#Mask out background gene sets
allgeneInfo_decay_gr<-mask_background(allgeneInfo_decay_gr)
#Mask out M-chromosome genes (except if chrM)
allgeneInfo_decay_gr<-mask_chromosome(allgeneInfo_decay_gr, chr_names=c('chrM'), feature_names_excluded=c("globalchromatin_M-chromosome"))
#Mask out Y-chromosome genes (except if chrY or chrM)
allgeneInfo_decay_gr<-mask_chromosome(allgeneInfo_decay_gr, chr_names=c('chrY'), feature_names_excluded=c("globalchromatin_Y-chromosome", "globalchromatin_M-chromosome"))
#Mask out X-chromosome genes for TAD-based features
allgeneInfo_decay_gr<-mask_chromosome(allgeneInfo_decay_gr, chr_names=c('chrX'), feature_names=colnames(elementMetadata(allgeneInfo_decay_gr))[grep("globalchromatin.*TAD", colnames(elementMetadata(allgeneInfo_decay_gr)))])

###########################
### Plot side-by-side comparison of chromatin gene set noise and mRNA stability bias
###########################

#Feature subset total variability (no repeats)
gst_chr_list_subset<-lapply(gst_chr_list, feature_subset)
plot_df<-do.call('rbind', gst_chr_list_subset)
plot_df$dataset<-rep(names(gst_chr_list_subset), times=as.numeric(sapply(gst_chr_list_subset, dim)[1,]))
gst_dec_all_list_subset<-lapply(gst_dec_all_list, feature_subset)
plot_df2<-do.call('rbind', gst_dec_all_list_subset)
plot_df2$dataset<-rep(names(gst_dec_all_list_subset), times=as.numeric(sapply(gst_dec_all_list_subset, dim)[1,]))
plot_df<-rbind(plot_df[,c("effect_size", "p_value", "FDR", "dataset")], plot_df2)
#Subset datasets
plot_df<-plot_df[plot_df$dataset %in% c('islam', 'grunS', 'grun2i', 'friedel', 'sharova', 'tippmann', 'schwanhausser'),]
#Dataset types (replicate datasets either of type 'expression noise' or 'mRNA stability')
plot_df$dataset_type<-'expression_noise'
plot_df$dataset_type[plot_df$dataset %in% c('friedel', 'sharova', 'tippmann', 'schwanhausser')]<-'mRNA_stability'
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
meds<-ddply(plot_df[plot_df$dataset_type=='expression_noise',], c("feat_category", "feat_name", "location"), summarise, effect_size_median = median(effect_size))
meds<-meds[order(meds$feat_category, abs(meds$effect_size_median-0.5), decreasing=T),]
meds<-meds[!duplicated(meds$feat_name),]
meds<-meds[order(meds$feat_category, meds$effect_size_median, decreasing=T),]
plot_df$index<-factor(plot_df$feat_name, levels=meds$feat_name)
#Plot
pdf('gene_set_variability_box_decsubset_ccc.pdf', width=10, height=7)
d <- ggplot(plot_df, aes(index, effect_size)) + 
  geom_boxplot(aes(color = location)) + geom_point(aes(color = location), size=1, position = position_dodge(width = 0.75)) +
  #scale_y_continuous(limits = c(0.35, 0.7)) +
  scale_colour_manual(values = plot_cols) + 
  theme_bw() + 
  geom_hline(yintercept = 0.5, size=0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d + facet_wrap(~dataset_type, ncol=2, nrow=1)
dev.off()

###########################
### Select mRNA stability dataset for downstream analyses
###########################

# #Select mRNA stability dataset (Sharova)
# decay_name<-'sharova'
# decay_dir<-paste('select', decay_name, sep='_')
# #Create output directory
# dir.create(decay_dir, showWarnings = FALSE)
# decay_select<-decay_sharova
# decay_select$mRNA_stability<-decay_sharova$bined.1
# #
#Select mRNA stability dataset (Friedel)
decay_name<-'friedel'
decay_dir<-paste('select', decay_name, sep='_')
#Create output directory
dir.create(decay_dir, showWarnings = FALSE)
decay_select<-decay_friedel
decay_select$mRNA_stability<-decay_friedel$Half.life..h.
# #
# #Select mRNA stability dataset (Schwanhausser)
# decay_name<-'schwanhausser'
# decay_dir<-paste('select', decay_name, sep='_')
# #Create output directory
# dir.create(decay_dir, showWarnings = FALSE)
# decay_select<-decay_schwanhausser
# decay_select$mRNA_stability<-decay_schwanhausser$mRNA.half.life.average..h.
# #
# #Select mRNA stability dataset (Tippmann)
# decay_name<-'tippmann'
# decay_dir<-paste('select', decay_name, sep='_')
# decay_select<-decay_tippmann
# decay_select$mRNA_stability<-decay_tippmann$HL_ES

###########################
### Plot mRNA half-life distributions for all chromosomes
###########################

#mRNA half-life distributions for all chromosomes
plot_df<-decay_select
plot_df$mRNA_stability<-plot_df$mRNA_stability
plot_df$chromosome<-as.character(seqnames(allgeneInfo_gr[rownames(decay_select),]))
plot_df$chromosome_type<-'autosome'
plot_df$chromosome_type[plot_df$chromosome=='chrX']<-'X-chromosome'
#Remove Y-chromosome
plot_df<-plot_df[plot_df$chromosome!='chrY',]
#Ordered factor
plot_df$chromosome<-factor(plot_df$chromosome, levels=paste('chr', c('X', 1:19), sep=''))
plot_cols<-c('grey', brewer.pal(5, 'Spectral')[5])
names(plot_cols)<-c('autosome', 'X-chromosome')
#Plot
pdf(paste(decay_dir, 'mRNA_hl_perchromosome_box.pdf', sep='/'), width=7, height=5)
d <- ggplot(plot_df, aes(chromosome, mRNA_stability)) + 
  geom_boxplot(aes(fill = chromosome_type), notch=T, outlier.shape = NA) + 
  scale_fill_manual(values = plot_cols) + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 21)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
d
dev.off()

###########################
### Plot mRNA half-life distributions of X chromosome compared to random sampling of matched size
###########################

#mRNA half-life distributions of X chromosome compared to random sampling of matched size
set.seed(1)
rand_repl<-1e4
rand_list<-list()
temp_genes<-rownames(decay_select)
temp_genes_X<-temp_genes[which(allgeneInfo_gr[temp_genes,]$"globalchromatin_X-chromosome"==1)]
temp_genes_noX<-temp_genes[which(allgeneInfo_gr[temp_genes,]$"globalchromatin_X-chromosome"==0)]
for(j in 1:rand_repl){
  rand_list[[j]]<-mean(decay_select[sample(temp_genes_noX, length(temp_genes_X), replace=F),]$mRNA_stability)
}
plot_df<-data.frame(mRNA_stability=unlist(rand_list))
#Plot
pdf(paste(decay_dir, 'mRNA_hl_chrX_hist_nullmodel.pdf', sep='/'), width=7, height=5)
d <- ggplot(plot_df, aes(mRNA_stability)) + 
  geom_histogram(colour='black', fill='grey') + 
  theme_bw() +
  geom_vline(xintercept = mean(decay_select[temp_genes_X,]$mRNA_stability), size=1, linetype=2, colour=brewer.pal(5, 'Spectral')[5])
d
dev.off()

###########################
### Weighted sampling of genes selected according to mRNA stability quantile
###########################

#Weighted sampling of genes selected according to mRNA stability quantile
set.seed(1)
decay_quantiles<-10
decay_quantile_replicates<-10
decay_quantile_sizes<-seq(200, 2000, 200)
temp_genes<-rownames(decay_select)
temp_quant<-quantile_fact(decay_select$mRNA_stability, 10)
dselect_wsample_list<-list()
for(j in 1:decay_quantiles){
  dselect_wsample_list[[j]]<-list()
  temp_quant_j<-decay_quantiles-abs(as.numeric(temp_quant)-j)
  for(k in 1:decay_quantile_replicates){
    for(l in 1:length(decay_quantile_sizes)){
      dselect_wsample_list[[j]][[(k-1)*10+l]]<-sample(temp_genes, decay_quantile_sizes[l], replace=F, prob=temp_quant_j)      
    }
  }
}

###########################
### Test biased decay of all decay-based gene sets
###########################

#Test biased decay of all decay-based gene sets
temp_decay<-decay_select$mRNA_stability
names(temp_decay)<-rownames(decay_select)
gsdt_dws_dselect<-genesets_biased(unlist(dselect_wsample_list, recursive=F), temp_decay)

###########################
### Test biased variance of decay quantiles
###########################

#Test biased variance of decay quantiles
gst_dws_dselect_list<-list() #biased expression variance of gene sets based on weighted sampling according to selected mRNA stability quantiles
temp_dirs<-list.dirs('..', recursive=F)
temp_dirs<-temp_dirs[grep('pagoda_adjvariance_bias', temp_dirs)]
#Retain cell cycle corrected data
temp_dirs<-temp_dirs[grep('ccc', temp_dirs)]
# temp_dirs<-temp_dirs[-grep("covariancecontrol|graf|klein|zeisel", temp_dirs)]
temp_file<-paste(temp_dirs, list.files(temp_dirs, pattern='ccc.RData'), sep='/')
names(temp_file)<-sapply(strsplit(temp_dirs, '_'), '[', 4)
temp_file<-temp_file[names(temp_file) %in% c('islam', 'grunS', 'grun2i')]
for(i in 1:length(temp_file)){
  print(names(temp_file)[i])
  load(temp_file[[i]])
  #Rank-based gene expression variation measure
  genevar_rank<-binx_nearesty(varinfo$avmodes, genevar[names(varinfo$avmodes)], bin_size=101)
  #Remove NAs
  genevar_rank<-genevar_rank[!is.na(genevar_rank)]
  #Test biased variance of all decay-based gene sets
  gst_dws_dselect_list[[names(temp_file)[i]]]<-genesets_biased(unlist(dselect_wsample_list, recursive=F), genevar_rank[names(genevar_rank) %in% rownames(decay_select)])
}

###########################
### Plot relationship between mRNA stability and gene expression noise (random and actual chromatin gene sets)
###########################

#Plot relationship between stability and gene expression noise (random)
plot_df1<-apply(do.call('cbind', lapply(gst_dws_dselect_list[c('islam', 'grunS', 'grun2i')], '[', 'effect_size')), 1, mean)
plot_df<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=plot_df1)
cors <- round(cor(plot_df$expr_noise, plot_df$mRNA_stability, method='spearman', use='pairwise.complete.obs'), 2)
#Relationship between stability and gene expression noise (empirical)
gst_chr_list_subset<-lapply(gst_chr_list, feature_subset)
temp_noise<-apply(do.call('cbind', lapply(gst_chr_list_subset[c('islam', 'grunS', 'grun2i')], '[', 'effect_size')), 1, mean)
temp_noise_se<-apply(do.call('cbind', lapply(gst_chr_list_subset[c('islam', 'grunS', 'grun2i')], '[', 'effect_size')), 1, sd)/sqrt(3)
temp_noise_sig<-as.numeric(apply(do.call('cbind', lapply(gst_chr_list_subset[c('islam', 'grunS', 'grun2i')], '[', 'FDR'))<0.1, 1, sum)==3)
gst_dec_all_list_subset<-lapply(gst_dec_all_list, feature_subset)
temp_stability<-gst_dec_all_list_subset[[decay_name]]$effect_size
temp_stability_sig<-as.numeric(gst_dec_all_list_subset[[decay_name]]$FDR<0.1)
plot_df_chr<-data.frame(mRNA_stability=temp_stability, mRNA_stability_sig=temp_stability_sig, expr_noise=temp_noise, expr_noise_sig=temp_noise_sig, expr_noise_min=temp_noise-temp_noise_se, expr_noise_max=temp_noise+temp_noise_se)
limits <- aes(ymax = expr_noise + temp_noise_se, ymin=expr_noise - temp_noise_se)
#Feature category
plot_df_chr$feat_category<-extract_feature_category(rownames(plot_df_chr))
#Feature name
plot_df_chr$feat_name<-extract_feature_name(rownames(plot_df_chr))
#Bubble colour
plot_df_chr$location<-'1 TSS'
plot_df_chr$location[grep('corePromInfo', rownames(plot_df_chr))]<-'2 core promoter'
plot_df_chr$location[grep('promInfo', rownames(plot_df_chr))]<-'3 promoter'
plot_df_chr$location[grep('allgeneInfo', rownames(plot_df_chr))]<-'4 whole gene body'
plot_df_chr$location<-as.character(plot_df_chr$location)
plot_cols<-c(brewer.pal(5, 'Spectral')[-3])
names(plot_cols)<-c('1 TSS', '2 core promoter', '3 promoter', '4 whole gene body')
#Plot
pdf(paste(decay_dir, 'mRNA_hl_vs_noise_genesets_nullmodel.pdf', sep='/'), width=7, height=5)
d <- ggplot(plot_df, aes(mRNA_stability, expr_noise)) +
  geom_point(colour='grey') +
  theme_bw() +
  geom_smooth(method="lm", colour='black', size=1, se=T) +
  annotate("text", label = paste('r =', cors), x = 0.65, y = 0.65) +
  geom_pointrange(data=plot_df_chr, aes(x = mRNA_stability, ymin = expr_noise_min, ymax = expr_noise_max, y = expr_noise, color = location), size=1) +
  geom_text(data=plot_df_chr, aes(label=feat_name), size=2) +
  scale_colour_manual(values = plot_cols) + 
  geom_hline(yintercept = 0.5, size=0.5, linetype=2) + 
  geom_vline(xintercept = 0.5, size=0.5, linetype=2)
d
# d + facet_wrap(~dataset, nrow=1, ncol=3)
dev.off()
#Only retain significant noise associations or X-chromosome
plot_df_chr_sig<-plot_df_chr[which((plot_df_chr$expr_noise_sig==1 & !is.na(plot_df_chr$mRNA_stability_sig)) | plot_df_chr$feat_name=='X-chromosome'),]
#Plot only significant chromatin features
pdf(paste(decay_dir, 'mRNA_hl_vs_noise_siggenesets_nullmodel.pdf', sep='/'), width=7, height=5)
d <- ggplot(plot_df, aes(mRNA_stability, expr_noise)) +
  geom_point(colour='grey') +
  theme_bw() +
  geom_smooth(method="lm", colour='black', size=1, se=T) +
  annotate("text", label = paste('r =', cors), x = 0.65, y = 0.65) +
  geom_pointrange(data=plot_df_chr_sig, aes(x = mRNA_stability, ymin = expr_noise_min, ymax = expr_noise_max, y = expr_noise, color = location), size=1) +
  geom_text(data=plot_df_chr_sig, aes(label=feat_name), size=2) +
  scale_colour_manual(values = plot_cols) + 
  geom_hline(yintercept = 0.5, size=0.5, linetype=2) + 
  geom_vline(xintercept = 0.5, size=0.5, linetype=2)
d
# d + facet_wrap(~dataset, nrow=1, ncol=3)
dev.off()
#Plot only global chromatin features
#Only retain global chromatin features
plot_df_chr_global<-plot_df_chr[plot_df_chr$feat_category=='Global chromatin',]
pdf(paste(decay_dir, 'mRNA_hl_vs_noise_globalgenesets_nullmodel.pdf', sep='/'), width=7, height=5)
d <- ggplot(plot_df, aes(mRNA_stability, expr_noise)) +
  geom_point(colour='grey') +
  theme_bw() +
  geom_smooth(method="lm", colour='black', size=1, se=T) +
  annotate("text", label = paste('r =', cors), x = 0.65, y = 0.65) +
  geom_pointrange(data=plot_df_chr_global, aes(x = mRNA_stability, ymin = expr_noise_min, ymax = expr_noise_max, y = expr_noise, color = location), size=1) +
  geom_text(data=plot_df_chr_global, aes(label=feat_name), size=2) +
  scale_colour_manual(values = plot_cols) + 
  geom_hline(yintercept = 0.5, size=0.5, linetype=2) + 
  geom_vline(xintercept = 0.5, size=0.5, linetype=2)
d
# d + facet_wrap(~dataset, nrow=1, ncol=3)
dev.off()

#Empirical P-values of three UMI-based (residual to lm >= that under null model?)
#Islam
plot_df_islam<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_islam)
(sum(emp_lm$residuals>=(gst_chr_list_subset[['islam']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
                          predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability']))))+1)/1001
#GrunS
plot_df_grunS<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_grunS)
(sum(emp_lm$residuals>=(gst_chr_list_subset[['grunS']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
                          predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability']))))+1)/1001
#Grun2i
plot_df_grun2i<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_grun2i)
(sum(emp_lm$residuals>=(gst_chr_list_subset[['grun2i']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
                          predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability']))))+1)/1001

#Noise effect size = Delta AUC (observed - expected under null model?)
#Islam
plot_df_islam<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_islam)
(gst_chr_list_subset[['islam']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
    predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability'])))
#GrunS
plot_df_grunS<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_grunS)
(gst_chr_list_subset[['grunS']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
    predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability'])))
#Grun2i
plot_df_grun2i<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_grun2i)
(gst_chr_list_subset[['grun2i']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
    predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability'])))

###########################
### Plot relationship between mRNA stability and gene expression noise (random and actual chromatin gene sets) => separately for each dataset
###########################

#Relationship between stability and gene expression noise (random) => separately for each dataset
for(focus_dataset in c('islam', 'grunS', 'grun2i')){
  plot_df<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[[focus_dataset]]$effect_size)
  cors <- round(cor(plot_df$expr_noise, plot_df$mRNA_stability, method='spearman', use='pairwise.complete.obs'), 2)
  #Relationship between stability and gene expression noise (empirical)
  gst_chr_list_subset<-lapply(gst_chr_list, feature_subset)
  temp_noise_sig<-as.numeric(apply(do.call('cbind', lapply(gst_chr_list_subset[c('islam', 'grunS', 'grun2i')], '[', 'FDR'))<0.1, 1, sum)==3)
  gst_dec_all_list_subset<-lapply(gst_dec_all_list, feature_subset)
  temp_stability<-gst_dec_all_list_subset[[decay_name]]$effect_size
  temp_stability_sig<-as.numeric(gst_dec_all_list_subset[[decay_name]]$FDR<0.1)
  plot_df_chr<-data.frame(mRNA_stability=temp_stability, 
                          mRNA_stability_sig=temp_stability_sig, 
                          expr_noise=gst_chr_list_subset[[focus_dataset]][,'effect_size'], 
                          expr_noise_sig=temp_noise_sig, row.names=rownames(gst_chr_list_subset[[focus_dataset]]))
  #Feature category
  plot_df_chr$feat_category<-extract_feature_category(rownames(plot_df_chr))
  #Feature name
  plot_df_chr$feat_name<-extract_feature_name(rownames(plot_df_chr))
  #Bubble colour
  plot_df_chr$location<-'1 TSS'
  plot_df_chr$location[grep('corePromInfo', rownames(plot_df_chr))]<-'2 core promoter'
  plot_df_chr$location[grep('promInfo', rownames(plot_df_chr))]<-'3 promoter'
  plot_df_chr$location[grep('allgeneInfo', rownames(plot_df_chr))]<-'4 whole gene body'
  plot_df_chr$location<-as.character(plot_df_chr$location)
  plot_cols<-c(brewer.pal(5, 'Spectral')[-3])
  names(plot_cols)<-c('1 TSS', '2 core promoter', '3 promoter', '4 whole gene body')
  #Plot
  d <- ggplot(plot_df, aes(mRNA_stability, expr_noise)) +
    geom_point(colour='grey') +
    theme_bw() +
    geom_smooth(method="lm", colour='black', size=1, se=F) +
    annotate("text", label = paste('r =', cors), x = 0.65, y = 0.65) +
    geom_point(data=plot_df_chr, aes(color = location), size=5) +
    geom_text(data=plot_df_chr, aes(label=feat_name), size=2) +
    scale_colour_manual(values = plot_cols) + 
    geom_hline(yintercept = 0.5, size=0.5, linetype=2) + 
    geom_vline(xintercept = 0.5, size=0.5, linetype=2)
  ggsave(paste(decay_dir, paste('mRNA_hl_vs_noise_genesets_nullmodel_', focus_dataset, '.pdf', sep=''), sep='/'), width=7, height=5)
  #Only retain significant noise associations or X-chromosome
  plot_df_chr_sig<-plot_df_chr[which((plot_df_chr$expr_noise_sig==1 & !is.na(plot_df_chr$mRNA_stability_sig)) | plot_df_chr$feat_name=='X-chromosome'),]
  #Plot only significant chromatin features
  pdf(paste(decay_dir, paste('mRNA_hl_vs_noise_siggenesets_nullmodel_', focus_dataset, '.pdf', sep=''), sep='/'), width=7, height=5)
  d <- ggplot(plot_df, aes(mRNA_stability, expr_noise)) +
    geom_point(colour='grey') +
    theme_bw() +
    geom_smooth(method="lm", colour='black', size=1, se=F) +
    annotate("text", label = paste('r =', cors), x = 0.65, y = 0.65) +
    geom_point(data=plot_df_chr_sig, aes(color = location), size=5) +
    geom_text(data=plot_df_chr_sig, aes(label=feat_name), size=2) +
    scale_colour_manual(values = plot_cols) + 
    geom_hline(yintercept = 0.5, size=0.5, linetype=2) + 
    geom_vline(xintercept = 0.5, size=0.5, linetype=2)
  ggsave(paste(decay_dir, paste('mRNA_hl_vs_noise_siggenesets_nullmodel_', focus_dataset, '.pdf', sep=''), sep='/'), width=7, height=5)
  #Plot only global chromatin features
  #Only retain global chromatin features
  plot_df_chr_global<-plot_df_chr[plot_df_chr$feat_category=='Global chromatin',]
  pdf(paste(decay_dir, paste('mRNA_hl_vs_noise_globalgenesets_nullmodel_', focus_dataset, '.pdf', sep=''), sep='/'), width=7, height=5)
  d <- ggplot(plot_df, aes(mRNA_stability, expr_noise)) +
    geom_point(colour='grey') +
    theme_bw() +
    geom_smooth(method="lm", colour='black', size=1, se=F) +
    annotate("text", label = paste('r =', cors), x = 0.65, y = 0.65) +
    geom_point(data=plot_df_chr_global, aes(color = location), size=5) +
    geom_text(data=plot_df_chr_global, aes(label=feat_name), size=2) +
    scale_colour_manual(values = plot_cols) + 
    geom_hline(yintercept = 0.5, size=0.5, linetype=2) + 
    geom_vline(xintercept = 0.5, size=0.5, linetype=2)
  ggsave(paste(decay_dir, paste('mRNA_hl_vs_noise_globalgenesets_nullmodel_', focus_dataset, '.pdf', sep=''), sep='/'), width=7, height=5)
}

#Empirical P-values of three UMI-based (residual to lm >= that under null model?)
#Islam
plot_df_islam<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_islam)
(sum(emp_lm$residuals>=(gst_chr_list_subset[['islam']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
                          predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability']))))+1)/1001
#GrunS
plot_df_grunS<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_grunS)
(sum(emp_lm$residuals>=(gst_chr_list_subset[['grunS']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
                          predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability']))))+1)/1001
#Grun2i
plot_df_grun2i<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_grun2i)
(sum(emp_lm$residuals>=(gst_chr_list_subset[['grun2i']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
                          predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability']))))+1)/1001

#Noise effect size = Delta AUC (observed - expected under null model?)
#Islam
plot_df_islam<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_islam)
(gst_chr_list_subset[['islam']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
    predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability'])))
#GrunS
plot_df_grunS<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_grunS)
(gst_chr_list_subset[['grunS']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
    predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability'])))
#Grun2i
plot_df_grun2i<-data.frame(mRNA_stability=gsdt_dws_dselect$effect_size, expr_noise=gst_dws_dselect_list[['islam']]$effect_size)
emp_lm<-lm(expr_noise~mRNA_stability, data=plot_df_grun2i)
(gst_chr_list_subset[['grun2i']]['allgeneInfo.globalchromatin_X-chromosome',]$effect_size-
    predict(emp_lm, data.frame(mRNA_stability=plot_df_chr['allgeneInfo.globalchromatin_X-chromosome','mRNA_stability'])))

