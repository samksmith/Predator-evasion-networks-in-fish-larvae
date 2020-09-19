# Network-Inference-Annotations
Samantha K. Smith, UT Austin, samksmith@utexas.edu

This method uses the output of gene expression analyses to create contextual and functional annotation of genes within relevant networks. We used this method on gene expression data from red drum larvae. This project contains scripts for generating annotations as well as the red drum behavioral data and associated annotation files. Red drum reads used in this study will be available at the Short Read Archive (NCBI SRA) soon.

# Input
NI_gene_annotation.R requires 2 files
1. Module csv files. Genes within the module have a listed module membership value (correlation between gene and module eigengene) for genes within the module and a "0" for all other genes.

2. Module-trait table. Lists significant correlations between modules and measured traits (module, trait, correlation, p-value).

# Output
A two-column csv file that contains gene ID and network inference annotation. Annotations are formatted as follows: 

   gene ID | module | kME value | Associations: (ME-trait correlation) trait1; (ME-trait correlation) trait2 

# Files
NI_gene_annotation.R - R script for creating network inference gene annotations

ge_pipeline_reddrum.R - Gene expression analysis pipeline used on red drum dataset. Majority of the pipeline was originally created by Mikhail Matz (https://github.com/z0on) and modified by Sam Smith.

Reddrum_annotations.csv - two column file with gene ID and network inference annotations for red drum larvae

reddrum_NIannotation_Eggnog.csv - three column file with gene ID, network inference annotations, and EggNOG-derived annotations

module-trait_NoPCs.csv - Lists significant correlations between modules and measured traits for red drum dataset (module, trait, correlation, p-value).

trait_to_name.csv - associates trait names with human readable trait descriptions

color_allrem.csv - module csv files for the red drum dataset. These modules were all significantly associated (p<=0.01) with at least one measured trait.
