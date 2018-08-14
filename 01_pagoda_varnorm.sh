#!/bin/bash

#Assumes current working directory is repository base path ("chromatin_noise_paper")

#Grun2i
mkdir pagoda_varnorm_grun2i
./pagoda_varnorm.R -i ./preQC_grun2i/preQC_grun2i.RData --maxTranscriptLengthFile ./misc_data_files/SC_preprocessing_transcriptmaxlength.RData --numCores 4 --minDetected 2 --minCountThreshold 1 --maxAdjVar 5 -o pagoda_varnorm_grun2i -p pagoda_varnorm_grun2i

#Grun2i pool-and-split
mkdir pagoda_varnorm_grun2iPandS
./pagoda_varnorm.R -i ./preQC_grun2iPandS/preQC_grun2iPandS.RData --maxTranscriptLengthFile ./misc_data_files/SC_preprocessing_transcriptmaxlength.RData --numCores 4 --minDetected 2 --minCountThreshold 1 --maxAdjVar 5 -o pagoda_varnorm_grun2iPandS -p pagoda_varnorm_grun2iPandS

#GrunS
mkdir pagoda_varnorm_grunS
./pagoda_varnorm.R -i ./preQC_grunS/preQC_grunS.RData --maxTranscriptLengthFile ./misc_data_files/SC_preprocessing_transcriptmaxlength.RData --numCores 4 --minDetected 2 --minCountThreshold 1 --maxAdjVar 5 -o pagoda_varnorm_grunS -p pagoda_varnorm_grunS

#GrunS pool-and-split
mkdir pagoda_varnorm_grunSPandS
./pagoda_varnorm.R -i ./preQC_grunSPandS/preQC_grunSPandS.RData --maxTranscriptLengthFile ./misc_data_files/SC_preprocessing_transcriptmaxlength.RData --numCores 4 --minDetected 2 --minCountThreshold 1 --maxAdjVar 5 -o pagoda_varnorm_grunSPandS -p pagoda_varnorm_grunSPandS

#Islam
mkdir pagoda_varnorm_islam
./pagoda_varnorm.R -i ./preQC_islam/preQC_islam.RData --maxTranscriptLengthFile ./misc_data_files/SC_preprocessing_transcriptmaxlength.RData --numCores 4 --minDetected 2 --minCountThreshold 1 --maxAdjVar 5 -o pagoda_varnorm_islam -p pagoda_varnorm_islam

#Teichmann2i
mkdir pagoda_varnorm_teichmann2i
./pagoda_varnorm.R -i ./preQC_teichmann2i/preQC_teichmann2i.RData --maxTranscriptLengthFile ./misc_data_files/SC_preprocessing_transcriptmaxlength.RData --numCores 4 --minDetected 2 --minCountThreshold 2 --maxAdjVar 5 -o pagoda_varnorm_teichmann2i -p pagoda_varnorm_teichmann2i

#Teichmanna2i
mkdir pagoda_varnorm_teichmanna2i
./pagoda_varnorm.R -i ./preQC_teichmanna2i/preQC_teichmanna2i.RData --maxTranscriptLengthFile ./misc_data_files/SC_preprocessing_transcriptmaxlength.RData --numCores 4 --minDetected 2 --minCountThreshold 2 --maxAdjVar 5 -o pagoda_varnorm_teichmanna2i -p pagoda_varnorm_teichmanna2i

#TeichmannS
mkdir pagoda_varnorm_teichmannS
./pagoda_varnorm.R -i ./preQC_teichmannS/preQC_teichmannS.RData --maxTranscriptLengthFile ./misc_data_files/SC_preprocessing_transcriptmaxlength.RData --numCores 4 --minDetected 2 --minCountThreshold 2 --maxAdjVar 5 -o pagoda_varnorm_teichmannS -p pagoda_varnorm_teichmannS
