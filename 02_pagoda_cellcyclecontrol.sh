#!/bin/bash

#Assumes current working directory is repository base path ("chromatin_noise_paper")

#Grun2i
mkdir pagoda_cellcyclecontrol_grun2i
./pagoda_cellcyclecontrol.R -i ./pagoda_varnorm_grun2i/pagoda_varnorm_grun2i.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_cellcyclecontrol_grun2i -p pagoda_cellcyclecontrol_grun2i

#Grun2i pool-and-split
mkdir pagoda_cellcyclecontrol_grun2iPandS
./pagoda_cellcyclecontrol.R -i ./pagoda_varnorm_grun2iPandS/pagoda_varnorm_grun2iPandS.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_cellcyclecontrol_grun2iPandS -p pagoda_cellcyclecontrol_grun2iPandS

#GrunS
mkdir pagoda_cellcyclecontrol_grunS
./pagoda_cellcyclecontrol.R -i ./pagoda_varnorm_grunS/pagoda_varnorm_grunS.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_cellcyclecontrol_grunS -p pagoda_cellcyclecontrol_grunS

#GrunS pool-and-split
mkdir pagoda_cellcyclecontrol_grunSPandS
./pagoda_cellcyclecontrol.R -i ./pagoda_varnorm_grunSPandS/pagoda_varnorm_grunSPandS.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_cellcyclecontrol_grunSPandS -p pagoda_cellcyclecontrol_grunSPandS

#Islam
mkdir pagoda_cellcyclecontrol_islam
./pagoda_cellcyclecontrol.R -i ./pagoda_varnorm_islam/pagoda_varnorm_islam.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_cellcyclecontrol_islam -p pagoda_cellcyclecontrol_islam

#Teichmann2i
mkdir pagoda_cellcyclecontrol_teichmann2i
./pagoda_cellcyclecontrol.R -i ./pagoda_varnorm_teichmann2i/pagoda_varnorm_teichmann2i.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_cellcyclecontrol_teichmann2i -p pagoda_cellcyclecontrol_teichmann2i

#Teichmanna2i
mkdir pagoda_cellcyclecontrol_teichmanna2i
./pagoda_cellcyclecontrol.R -i ./pagoda_varnorm_teichmanna2i/pagoda_varnorm_teichmanna2i.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_cellcyclecontrol_teichmanna2i -p pagoda_cellcyclecontrol_teichmanna2i

#TeichmannS
mkdir pagoda_cellcyclecontrol_teichmannS
./pagoda_cellcyclecontrol.R -i ./pagoda_varnorm_teichmannS/pagoda_varnorm_teichmannS.RData --GOFile ./misc_data_files/GOdata_3ontologies.RData --numCores 4 --nRandomizations 50 --maxPathwaySize 2000 -o pagoda_cellcyclecontrol_teichmannS -p pagoda_cellcyclecontrol_teichmannS
