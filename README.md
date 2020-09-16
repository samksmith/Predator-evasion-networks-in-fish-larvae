# Network-Inference-Annotations
Samantha K. Smith, UT Austin, samksmith@utexas.edu

This method uses the output of gene expression analyses to create contextual and functional annotation of genes within relevant networks. We used this method on gene expression data from red drum larvae. This project contains scripts for generating annotations as well as the red drum behavioral data and associated annotation files. Red drum reads used in this study will be available at the Short Read Archive (NCBI SRA: https://www-ncbi-nlm-nih-gov.ezproxy.lib.utexas.edu/sra) soon.

# Input
NI_gene_annotation.R requires 2 files
1. Module csv files. Genes within the module have a listed module membership value (correlation between gene and module eigengene) for genes within the module and a "0" for all other genes.

2. Module-trait table. Lists significant correlations between modules and measured traits (module, trait, correlation, p-value).

# Output
A two-column csv file that contains gene ID and network inference annotation. Annotations are formatted as follows: 

   gene ID | module | kME value | Associations: (ME-trait correlation) trait1; (ME-trait 288correlation) trait2  
