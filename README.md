# Code for the Paraphelidium transcriptome (Torruella et al., Nat Comm Biol 2018)

This repository contains scripts for the *Paraphelidium tribonemae* transcriptome paper (Torruella et al., Nat Comm Biol 2018). It contains the basic R code to analyse the profile of presence/absence of COG functional annotations across various eukaryotic genomes, using the following methods:

* Principal Coordinate Analysis
* Clustering of species (based on Pearson's correlation coefficient and Ward clustering) 

For a detailed overview of the methods, see Methods in the paper [link].

## Required libraries

* Base *R*, includinng *stats* library
* *gplots*, for heatmap plots
* *ape*, for PCoA  analysis

## Input

This script can be manually edit to use various matrices as input. They should have the following format:
* Each column is a species
* Each row is a functional annotation group, e.g. clusters of orthologous groups (COG)
* The matrix will contain the presence/absence pattern of a given COG (encoded as 1/0).

The script can also be used to convert count matrices (e.g. number of genes with a certain COG) into presence/abscence matrices (1/0).

## Output

This code can be used to produce figures analogous to Figure 3 in the paper [link].
