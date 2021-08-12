# Predator-evasion-networks-in-fish-larvae
Samantha K. Smith, UT Austin, samksmith@utexas.edu

This is the repository of data and scripts for a forthcoming paper analyzing RNA sequencing data from red drum larvae that underwent startle assays. Red drum reads used in this study will be available at the Short Read Archive (NCBI SRA) soon.

# Input
Network_gene_annotation.R requires 2 types of files.
1. Module csv files. Genes within the module have a listed module membership value (correlation between gene and module eigengene) for genes within the module and a "0" for all other genes.

2. Module-trait table. Lists significant correlations between modules and measured traits (module, trait, correlation, p-value).

# Output
A two-column csv file that contains gene ID and annotation. Annotations are formatted as follows: 

   gene ID | module | kME value | Associations: (ME-trait correlation) trait1; (ME-trait correlation) trait2 

# Files

allcounts_redDrum.txt - raw counts data for red drum dataset.

color_allrem.csv - module csv files for the red drum dataset. These modules were all significantly associated (p<=0.01) with at least one measured trait.

EggNOG2GO.txt - Pipeline for creating eggNOG gene GO annotations.

fuimanfish_metadata.csv - One of three datasets containing measured traits and information about red drum larvae.

ge_pipeline_reddrum.R - Gene expression analysis pipeline used on red drum dataset. Majority of the pipeline was originally created by Mikhail Matz (https://github.com/z0on) and modified by Sam Smith.

genetic_relatedness.R - Pipeline for detecting relatedness between red drum larvae. Majority of the pipeline was originally created by Mikhail Matz (https://github.com/z0on) and modified by Sam Smith.

mito_traits.covMat - covariance matrix based on single read resampling for mitochondrial genes in the red drum dataset.

mito_traits.ibsMat - identity-by-state matrix for mitochondrial genes in the red drum dataset.

module-trait_NoPCs.csv - Lists significant correlations between modules and measured traits for red drum dataset (module, trait, correlation, p-value).

mom_id_fuimanfish.txt - Pipeline for making covariance and identity-by-state matrices for mitochondrial genes and all genes.

myresults.covMat - covariance matrix based on single read resampling for all genes in the red drum dataset.

myresults.ibsMat - identity-by-state matrix for all genes in the red drum dataset.

performance_reorg_042120.csv - One of three datasets containing measured traits from red drum larvae dataset.

module-trait_NoPCs.csv - List of significant correlations between modules and traits for the red drum dataset.

traits_and_time.csv - One of three datasets containing measured traits from red drum. Includes timing between behavioral assay and tissue collection.

trait_to_name.csv - associates trait names with human readable trait descriptions
