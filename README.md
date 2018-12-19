# Code for the Paraphelidium transcriptome (Torruella et al., Nat Comm Biol 2018)

This repository contains scripts for the *Paraphelidium tribonemae* transcriptome paper ([Torruella et al., Nat Comm Biol 2018](https://www.nature.com/articles/s42003-018-0235-z)). It contains the basic R code to analyse the profile of presence/absence of genes with certain functional annotations across various eukaryotic genomes, using the following methods:

* Principal Coordinate Analysis
* Clustering of species (based on Pearson's correlation coefficient and Ward clustering) 

In the paper, we used this script to analyse the profiles of clusters of orthologous groups (COGs) linked to primary metabolism (COGs: codes: C, E, F, G, H, I, P and Q), and KEGG orthologs (KOGs). For a detailed overview of the functional annotations of genes with COG and KEGG terms, see [Methods in the paper](https://www.nature.com/articles/s42003-018-0235-z#Sec8).

## Required libraries

* Base *R*, includinng *stats* library
* *gplots*, for heatmap plots
* *ape*, for PCoA  analysis
* *viridis* color palette
* *reshape2* and *tidyr* 

## Input

This script can be manually edited to use various matrices as input. They should have the following format:
* Each column is a species
* Each row is a functional annotation group, e.g. clusters of orthologous groups (KOG)
* The matrix will can contain either i) counts of a certain gene family (COG, KOG, etc.) in each species; or ii) presence/absence profile (encoded as 1/0). If counts are used, the script will automatically convert them to 0/1

The input file provided, ``primary_transposed.txt``, contains counts of KOGs across 41 eukaryotes (as in Figure 3 from the manuscript).

## Output

This code can be used to produce figures analogous to [Figure 3 in the paper](https://www.nature.com/articles/s42003-018-0235-z#Fig3).

## Cite this

Torruella G, Grau-Bové X, Moreira D, Karpov SA, Burns JA, Sebé-Pedrós A, Völcker E, López-García P. 2018. Global transcriptome analysis of the aphelid Paraphelidium tribonemae supports the phagotrophic origin of fungi. Commun Biol 1: 231. http://www.nature.com/articles/s42003-018-0235-z.
