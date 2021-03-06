# Brassica oleracea/rapa/napus code
Code for the Brassica oleracea/rapa/napus genomic comparison.

Most of these folders are RStudio projects that can be opened in RStudio.

## 0-Assembly_Annotation

This folder collects small scripts and documentation as to how I ran the assembly and annotation steps for the pangenome.

## 1-FilterPAVTables

This contains the PAV tables and some code to filter them to be species-specific and to remove 0 individual genes. 

## 2-PanPanComparison

This folder contains scripts to compare gene content between the pangenomes

## 3-Modeling

R-code for pangenome size modeling

## 4-VennDiagrams

R-code for the Venn-diagram plots 

## 5-Cluster_Napus_A_and_C_PAV

R-code for the PCA plots and animations 

## 6-FPSc

Code to compare the B. rapa FPSc individuals with the rest lives in this folder 

## 7-CompareRGenes and R_Gene_Analysis

R-code for the stacked barplots for the R-genes

## 8-STRINGS_Visualisation

Code to visualise differences in networked genes

## 9-PAVModeling

PAV modeling scripts 

## 10-GOEnrichment

Code and RStudio project that checks which GO terms are enriched

## 11-GenesUniqueToEachPangenome

More PAV visualisation scripts

## 12-SupplementaryTables

Just some data with coverages per individual

## 13-Visualise_Core_A_C

R-code for the variable genes vs core genes

## 14-FindIncompatibleGenes

Code to find mutually incompatible genes

## 15-SHAP_replotting

R-code to replot SHAP plots with defaults and values I like, one line per chromosome, a nicer grid and so on

## 16-ChisquareTestVariableGenesNetworks

A nicer recode of older code that runs chisquare for a bunch of comparison tables around variable/core genes in protein-protein interaction networks.

## 17-makeBetterSupplTables

R-code that pulls out coverage and other stats from SRA for the individuals I used.

## Notebooks 
These contain the XGBoost and Shapley values work used for modeling gene presence/absence variation.

You can launch them on binder:
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/appliedbioinformatics/Brassica_oleracea_rapa_napus_code/master?filepath=notebooks)

Or navigate to the `notebooks` folder.
