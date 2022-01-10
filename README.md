# AltSplDevCancer
#This repository contains the code used to identify embryonic splicing events and assess their relevence to the cellular physiology in cancer patients. 

A breif description of the scripts is provided below:

The script DownloadTrascriptTPMdataFromXena.R was used to download the transcript level TPMs from UCSC xena browser for the relevant tissues and 
cancer types of interest.

The script KallistoPipelineArrayExpressForSlurm.sh was used to download and process the raw RNA-seq data for human tissues from GEO database and was run 
using NIH highthroughput computing environment called biowulf.

The script ScorePathwaysIndevelopment.R uses principal component analuysis of the temporal transcriptomic data of human organs to identify the kegg pathways 
which are preferetially active during the pre-natal stage. Then it annotates the alternatively spliced exons as embryonic positive (EP) and 
embryonic negative (EN) based on the their coorrelation with embryonic pathways.


The script FrequentOutlierSplicingEventsInCancerPatients.R was used to detect the splicing events which frequently increased or deacresed 
(i.e deviated from corresponding Gtex normals by 2 standard deviations in atleast 20% of the cancer patients).
