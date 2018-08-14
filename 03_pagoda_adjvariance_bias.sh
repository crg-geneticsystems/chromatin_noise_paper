#!/bin/bash

#Assumes current working directory is repository base path ("chromatin_noise_paper")

#Grun2i
mkdir pagoda_adjvariance_bias_grun2i_ccc
./pagoda_adjvariance_bias.R -i ./pagoda_cellcyclecontrol_grun2i/pagoda_cellcyclecontrol_grun2i.RData --geneSetFile ./misc_data_files/SC_preprocessing_genesets.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_adjvariance_bias_grun2i_ccc -p pagoda_adjvariance_bias_grun2i_ccc

#Grun2i pool-and-split
mkdir pagoda_adjvariance_bias_grun2iPandS_ccc
./pagoda_adjvariance_bias.R -i ./pagoda_cellcyclecontrol_grun2iPandS/pagoda_cellcyclecontrol_grun2iPandS.RData --geneSetFile ./misc_data_files/SC_preprocessing_genesets.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_adjvariance_bias_grun2iPandS_ccc -p pagoda_adjvariance_bias_grun2iPandS_ccc

#GrunS
mkdir pagoda_adjvariance_bias_grunS_ccc
./pagoda_adjvariance_bias.R -i ./pagoda_cellcyclecontrol_grunS/pagoda_cellcyclecontrol_grunS.RData --geneSetFile ./misc_data_files/SC_preprocessing_genesets.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_adjvariance_bias_grunS_ccc -p pagoda_adjvariance_bias_grunS_cc

#GrunS pool-and-split
mkdir pagoda_adjvariance_bias_grunSPandS_ccc
./pagoda_adjvariance_bias.R -i ./pagoda_cellcyclecontrol_grunSPandS/pagoda_cellcyclecontrol_grunSPandS.RData --geneSetFile ./misc_data_files/SC_preprocessing_genesets.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_adjvariance_bias_grunSPandS_ccc -p pagoda_adjvariance_bias_grunSPandS_ccc

#Islam
mkdir pagoda_adjvariance_bias_islam_ccc
./pagoda_adjvariance_bias.R -i ./pagoda_cellcyclecontrol_islam/pagoda_cellcyclecontrol_islam.RData --geneSetFile ./misc_data_files/SC_preprocessing_genesets.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_adjvariance_bias_islam_ccc -p pagoda_adjvariance_bias_islam_ccc

#Teichmann2i
mkdir pagoda_adjvariance_bias_teichmann2i_ccc
./pagoda_adjvariance_bias.R -i .pagoda_cellcyclecontrol_teichmann2i/pagoda_cellcyclecontrol_teichmann2i.RData --geneSetFile ./misc_data_files/SC_preprocessing_genesets.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_adjvariance_bias_teichmann2i_ccc -p pagoda_adjvariance_bias_teichmann2i_ccc

#Teichmanna2i
mkdir pagoda_adjvariance_bias_teichmanna2i_ccc
./pagoda_adjvariance_bias.R -i ./pagoda_cellcyclecontrol_teichmanna2i/pagoda_cellcyclecontrol_teichmanna2i.RData --geneSetFile ./misc_data_files/SC_preprocessing_genesets.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_adjvariance_bias_teichmanna2i_ccc -p pagoda_adjvariance_bias_teichmanna2i_ccc

#TeichmannS
mkdir pagoda_adjvariance_bias_teichmannS_ccc
./pagoda_adjvariance_bias.R -i ./pagoda_cellcyclecontrol_teichmannS/pagoda_cellcyclecontrol_teichmannS.RData --geneSetFile ./misc_data_files/SC_preprocessing_genesets.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_adjvariance_bias_teichmannS_ccc -p pagoda_adjvariance_bias_teichmannS_ccc
