# Predator-evasion-networks-in-fish-larvae
Samantha K. Smith, UT Austin, samksmith@utexas.edu

This is the repository of data and scripts for a forthcoming paper analyzing RNA sequencing data from red drum larvae that underwent startle assays. Red drum reads used in this study will be available at the Short Read Archive (NCBI SRA) soon.

# Files

allcounts_redDrum.txt - raw counts data for red drum dataset.

color_allrem.csv - module csv files for the red drum dataset. These modules were all significantly associated (p<=0.01) with at least one measured trait and were enriched for GO terms (after removal of covariates and gpc1b from dataset).

EggNOG2GO.txt - Pipeline for creating eggNOG gene GO annotations.

fuimanfish_metadata.csv - One of three datasets containing measured traits and information about red drum larvae.

ge_pipeline_reddrum.R - Gene expression analysis pipeline used on red drum dataset. Pipeline was originally created by Mikhail Matz (https://github.com/z0on) and modified by Sam Smith.

genetic_relatedness.R - Pipeline for detecting relatedness between red drum larvae. Pipeline was originally created by Mikhail Matz (https://github.com/z0on) and modified by Sam Smith.

mito.covMat - covariance matrix based on single read resampling for mitochondrial genes in the red drum dataset.

mito.ibsMat - identity-by-state matrix for mitochondrial genes in the red drum dataset.

relatedness.txt - Pipeline for making covariance and identity-by-state matrices for mitochondrial genes and all genes.

allgenes.covMat - covariance matrix based on single read resampling for all genes in the red drum dataset.

allgenes.ibsMat - identity-by-state matrix for all genes in the red drum dataset.

performance_reorg_042120.csv - One of three datasets containing measured traits from red drum larvae dataset.

traits_and_time.csv - One of three datasets containing measured traits from red drum. Includes timing between behavioral assay and tissue collection.

trait_to_name.csv - associates trait names with human readable trait descriptions
