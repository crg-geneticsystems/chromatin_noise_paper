
#Determine desired quantiles of numeric input vector
#Returns factor vector of corresponding quantiles
quantile_fact<-function(in_vect, num_quant){
  temp_quant<-quantile(in_vect, probs=seq(0, 1, 1/num_quant))
  temp_fact<-rep(1, length(in_vect))
  for(i in 3:(num_quant+1)){
    temp_fact[in_vect<=temp_quant[i] & in_vect>temp_quant[i-1]]<-i-1
  }
  as.factor(temp_fact)
}

#Wrapper for Mann-Whitney U test
#Returns named vector with effect size (AUC) and p-value obtained from coin package
mann_whitney_U_wrapper<-function(vals1, vals2){
  if(length(vals1)==0 | length(vals2)==0){
    temp_test<-c(NA, NA)
  }else{
    g = factor(c(rep("GroupA", length(vals1)), rep("GroupB", length(vals2))))
    v = c(vals1, vals2)
    g_wt<-wilcox_test(v ~ g)
    temp_test<-c(wilcox.test(vals1, vals2)$statistic/(as.numeric(length(vals1))*as.numeric(length(vals2))),
                 pvalue(g_wt))
  }
  names(temp_test)<-c('effect_size', 'p_value')
  temp_test    
}

#Perform Mann-Whitney U test on numeric vectors in input lists (element indeces correspond)
#If background list not supplied, perform Mann-Whitney U test on each pair of numeric vectors in input list (first element in pair considered test vector)
#Returns data frame with effect size (AUC) and p-value obtained from coin package
mann_whitney_U_list<-function(test_list, background_list=NULL){
  if(is.null(background_list)){
    background_list<-test_list[seq(2, length(test_list), 2)]
    test_list<-test_list[seq(1, length(test_list), 2)]
  }
  temp_test<-as.data.frame(t(mapply(mann_whitney_U_wrapper, test_list, background_list)))
  temp_test$FDR<-p.adjust(temp_test$p_value)
  temp_test
}

#Ensembl gene annotations using data obtained from API
#Returns list of Ensembl genes with gene biotype (data frame) and corresponding GRanges object for chromosome subset of interest
ensemblAPI_genes_GRanges<-function(ens_file, chromosome_subset){
  temp<-read.table(ens_file, header=F, sep='\t', stringsAsFactors=F)
  colnames(temp)<-c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand', 'external_name', 'mgi_symbol', 'gene_biotype')
  allgenes<-temp
  #Reshuffle
  allgenes<-unique(allgenes[,c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand', 'gene_biotype', 'external_name')])
  allgenes<-allgenes[allgenes$chromosome_name %in% chromosome_subset,]
  rownames(allgenes)<-allgenes$ensembl_gene_id
  #Granges
  allgeneInfo_gr<-GRanges(seqnames=as(gsub('chrMT', 'chrM', paste('chr', allgenes$chromosome_name, sep='')), "Rle"), 
                          ranges=IRanges(allgenes$start_position, allgenes$end_position),
                          strand=as(allgenes$strand, "Rle"))
  names(allgeneInfo_gr)<-allgenes$ensembl_gene_id
  allgenes<-temp
  #Reshuffle
  allgenes<-unique(allgenes[,c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand', 'gene_biotype', 'external_name', 'mgi_symbol')])
  allgenes<-allgenes[allgenes$chromosome_name %in% chromosome_subset,]
  #Return
  list(allgenes=allgenes, allgeneInfo_gr=allgeneInfo_gr)
}

#Extract values from numeric input vector when index vector equals specified value
#Returns a numeric vector of values that satisfy the conditions specified by index_vect and index_val
index_vect<-function(num_vect, index_vect, index_val){
  num_vect[which(index_vect==index_val)]
}

#Test biased values for all gene sets (supplied as list)
#gene_sets_list either contains binary (double) vector matching num_vect or character vector matching names of num_vect
#Returns a data frame of results from the mann_whitney_U_list function
genesets_biased<-function(gene_sets_list, num_vect){
  if(typeof(gene_sets_list[[1]])=='double'){
    temp_values<-as.list(as.data.frame(matrix(num_vect, nrow=length(num_vect), ncol=length(gene_sets_list), byrow=F)))
    temp_in_list<-mapply(index_vect, temp_values, gene_sets_list, 1)
    temp_out_list<-mapply(index_vect, temp_values, gene_sets_list, 0)
    names(temp_in_list)<-names(gene_sets_list)  
    names(temp_out_list)<-names(gene_sets_list)
    mann_whitney_U_list(temp_in_list, temp_out_list)    
  }else if(typeof(gene_sets_list[[1]])=='character'){
    temp_in_list<-lapply(gene_sets_list, function(x){num_vect[names(num_vect) %in% x]})
    temp_out_list<-lapply(gene_sets_list, function(x){num_vect[!names(num_vect) %in% x]})
    mann_whitney_U_list(temp_in_list, temp_out_list)
  }
}

#Extract feature name
#Returns a feature name string
extract_feature_name<-function(char_vect){
  temp<-sapply(strsplit(char_vect, 'Info.'), '[', 2)
  temp[grep('sequence_core_promoter_', temp)]<-paste(' ', gsub('sequence_core_promoter_', '', temp[grep('sequence_core_promoter_', temp)]), sep='')
  temp[grep('sequence_repeat_', temp)]<-gsub('sequence_repeat_', '', temp[grep('sequence_repeat_', temp)])
  temp[grep('sequence_', temp)]<-gsub('sequence_', '', temp[grep('sequence_', temp)])
  temp[grep('localchromatin_peak_', temp)]<-gsub('localchromatin_peak_', '', temp[grep('localchromatin_peak_', temp)])
  temp[grep('localchromatin_state_', temp)]<-gsub('localchromatin_state_', '', temp[grep('localchromatin_state_', temp)])
  temp[grep('globalchromatin_', temp)]<-gsub('globalchromatin_', '', temp[grep('globalchromatin_', temp)])
  temp
}

#Extract feature category
#Returns a feature category string
extract_feature_category<-function(char_vect){
  temp<-''
  temp[grep('sequence_', char_vect)]<-'Sequence'
  temp[grep('sequence_repeat_', char_vect)]<-'Sequence (repeats)'
  temp[grep('localchromatin_peak_', char_vect)]<-'Local chromatin (peaks)'
  temp[grep('localchromatin_state_', char_vect)]<-'Local chromatin (state)'
  temp[grep('globalchromatin_', char_vect)]<-'Global chromatin'
  temp
}

#Replace GRanges metadata with NAs if a corresponding background feature exists (i.e. with "_background" suffix) and is 0
mask_background<-function(granges_object){
  temp_names<-colnames(elementMetadata(granges_object))
  temp_bgnames<-temp_names[grep('_background', temp_names)]
  temp_nbgnames<-temp_names[-grep('_background', temp_names)]
  for(i in temp_bgnames){
    elementMetadata(granges_object)[,gsub('_background', '', i)][elementMetadata(granges_object)[,i]==0]<-NA
  }
  elementMetadata(granges_object)<-elementMetadata(granges_object)[,temp_nbgnames]
  return(granges_object)
}

#Replace GRanges metadata with NAs if present on the specified chromosomes
mask_chromosome<-function(granges_object, chr_names, feature_names=NULL, feature_names_excluded=NULL){
  temp_names<-feature_names
  if(is.null(feature_names)){
    temp_names<-colnames(elementMetadata(granges_object))
  }
  if(!is.null(feature_names_excluded)){
    temp_names<-temp_names[!temp_names %in% feature_names_excluded]
  }
  for(i in temp_names){
    elementMetadata(granges_object)[,i][as.character(seqnames(granges_object)) %in% chr_names]<-NA
  }
  return(granges_object)
}

#Replace GRanges metadata with NAs if ranges less than the specified minimum length
mask_minlength<-function(granges_object, min_length, feature_names=NULL, feature_names_excluded=NULL){
  temp_names<-feature_names
  if(is.null(feature_names)){
    temp_names<-colnames(elementMetadata(granges_object))
  }
  if(!is.null(feature_names_excluded)){
    temp_names<-temp_names[!temp_names %in% feature_names_excluded]
  }
  for(i in temp_names){
    elementMetadata(granges_object)[,i][width(granges_object) < min_length]<-NA
  }
  return(granges_object)
}

#Subset gene set enrichment results to a custom set of features
feature_subset<-function(gs_tests_all, filter_pattern=NULL){
  if(!is.null(filter_pattern)){
    gs_tests_all<-gs_tests_all[grep(filter_pattern, gs_tests_all$name),]
  }
  gs_tests_all_subglobal<-gs_tests_all[c("allgeneInfo.globalchromatin_X-chromosome",
                                         "allgeneInfo.globalchromatin_Centromeric band",
                                         "allgeneInfo.globalchromatin_Telomeric band",
                                         "allgeneInfo.globalchromatin_cLAD",
                                         "allgeneInfo.globalchromatin_fLAD",
                                         "allgeneInfo.globalchromatin_Super-enhancer ES (target)"),]
  gs_tests_all_sublocal<-gs_tests_all[grep('promInfo.localchromatin_peak_H3k09me3|allgeneInfo.localchromatin_peak_H3k09me3', rownames(gs_tests_all)),]
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_H3k27me3', rownames(gs_tests_all)),])
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_H3k27ac', rownames(gs_tests_all)),])
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_H3k04me1|peak_H3k4me1', rownames(gs_tests_all)),])
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_H3k04me3|peak_H3k4me3', rownames(gs_tests_all)),])
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_H3k09ac', rownames(gs_tests_all)),])
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_H3k36me3', rownames(gs_tests_all)),])
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_P300', rownames(gs_tests_all)),])
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_Pol2', rownames(gs_tests_all)),])
  gs_tests_all_sublocal<-rbind(gs_tests_all_sublocal, gs_tests_all[grep('peak_Ctcf', rownames(gs_tests_all)),])
  gs_tests_all_subpromtype<-gs_tests_all[grep('CpG', rownames(gs_tests_all)),]
  gs_tests_all_subpromtype<-rbind(gs_tests_all_subpromtype, gs_tests_all[grep('Broad', rownames(gs_tests_all)),])
  gs_tests_all_subpromtype<-rbind(gs_tests_all_subpromtype, gs_tests_all[grep('CCAAT motif', rownames(gs_tests_all)),])
  gs_tests_all_subpromtype<-rbind(gs_tests_all_subpromtype, gs_tests_all[grep('GC motif', rownames(gs_tests_all)),])
  gs_tests_all_subpromtype<-rbind(gs_tests_all_subpromtype, gs_tests_all[grep('Initiator motif', rownames(gs_tests_all)),])
  gs_tests_all_subpromtype<-rbind(gs_tests_all_subpromtype, gs_tests_all[grep('Multiple initiation sites', rownames(gs_tests_all)),])
  gs_tests_all_subpromtype<-rbind(gs_tests_all_subpromtype, gs_tests_all[grep('Sharp TSS', rownames(gs_tests_all)),])
  gs_tests_all_subpromtype<-rbind(gs_tests_all_subpromtype, gs_tests_all[grep('TATAbox', rownames(gs_tests_all)),])
  #   gs_tests_all_subpromtype<-rbind(gs_tests_all_subpromtype, gs_tests_all[grep('Short', rownames(gs_tests_all)),])
  gs_tests_all_sub<-rbind(gs_tests_all_subglobal, gs_tests_all_sublocal, gs_tests_all_subpromtype)
  #   gs_tests_all_sub$FDR<-p.adjust(gs_tests_all_sub$p_value, method='BH')
  gs_tests_all_sub
}

#Rank-based normalisation
binx_nearesty<-function(var_x, var_y, bin_size=101){
  #Reorder
  temp_names<-names(var_y)
  temp_order<-1:length(var_x)
  temp_order<-temp_order[order(var_x)]
  var_y<-var_y[order(var_x)]
  var_x<-var_x[order(var_x)]
  #Indices of nearest n points (focus is first in list)
  var_xi<-lapply(as.list(1:length(var_x)), function(x){c(x, seq(x-(bin_size-1)/2, x+(bin_size-1)/2)[-((bin_size-1)/2+1)])})
  #Rank of focus within bin
  var_ybr<-sapply(lapply(lapply(var_xi, function(x){var_y[x[x>0]]}), rank), '[', 1)
  #Set incomplete bins to NA
  var_ybr[sapply(var_xi, min)<1 | sapply(var_xi, max)>length(var_xi)]<-NA
  #Reorder
  var_ybr<-var_ybr[order(temp_order)]
  #Rename
  names(var_ybr)<-temp_names
  #Return
  var_ybr
}

#Convert GRanges binary metadata columns indicating gene set membership to gene sets
binaryvect_to_geneset<-function(granges_object, name_prefix){
  temp_object<-list()
  for(i in colnames(elementMetadata(granges_object))){
    temp_object[[paste(name_prefix, 'yes', i, sep='_')]]<-na.omit(names(granges_object)[elementMetadata(granges_object)[,i]==1])
    temp_object[[paste(name_prefix, 'no', i, sep='_')]]<-na.omit(names(granges_object)[elementMetadata(granges_object)[,i]==0]) 
  }
  return(temp_object)
}

#Normalise features (columns in a data frame) prior to modelling
glm_prenorm<-function(x){
  for(i in 1:dim(x)[2]){
    if(length(unique(na.omit(x[,i])))>2 & !is.factor(x[,i])){
      x[,i]<-scale(x[,i])
    }
  }
  x
}

