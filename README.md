# Lepidoptera Orthologue Amino Acid Comparison

## Introduction
This repository is a collection of code and guides created by Arkan De Lomas for using orthofinder in bioinformatics. 
More information on Orthofinder can be found in its repository: https://github.com/davidemms/OrthoFinder

## Process
30+ lepidopteran species have their exomes analysed (using orthofinder) for gene sequence similarity to then be gathered into orthologues and non-orthologues. 
The data produced by Orthofinder contains a distinction between genes in each species that are orthologous and genes that are not (nonorthologous genes do not possess counterparts in other species within the sample).
amino acid content of orthologous and non-orthologous genes are then counted to produce frequency comparisons for each species to then deduce trends and niches within each species based on their climate, diet, and amino acid content of genes.

## Description of files in repository:
The main Orthogroup python file is Orthogroup_counter.py while for non-orthologous genes there are two python files as Orthofinder organizes the non-orthologous gene results differently from the orthologous gene results.

### Orthogroup_counter.py
This file produces a list of amino acid frequency in each othrologous gene for each species. 

### nonortholog_seqs.py
This file collates non-orth seq data from multiple columns into one colum and gathers frequency of amino acids

### non_orth_summary.py
Reads data produced by nonortholog_seqs.py ('nonorth_seq_and_aa.csv') and assigns them to species that contain those genes. 

