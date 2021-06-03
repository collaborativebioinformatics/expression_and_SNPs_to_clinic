# File: Combine.R
# Purpose: Integrate results from Variant Calling and Differential Expression
#          workflows into a unified tsv file for each sample.
#          
# Usage: 
# Author: Adam Faranda
# Created: June 3, 2021
###################################################################################


###################### Setup Environment and Import Tables ########################
library(dplyr)

### GATK Variants from DNA Seq

dna_seq_var

### GATK Variants from RNA Seq

rna_seq_var

### Differential Expression (Tumor vs Normal)

deg_table

######################### Generate Unified Identifier Table ##################

### Table consists of the Union of all distinct mutations observed 
### in the patient's Normal and Tumor samples
###
### Table Specification:
###     Variant HGVS expression
###     Hugo gene identifier
###     Variant sequence ontology consequence
###     ClinVar annotation
###     COSMIC annotation
###     Variant zygosity
###     Gene-level differential expression value (logFC, FDR, Abundances)
###     Source (DNA,RNA,Both)
###     Source identifier / Tissue origin, if relevant (Essential)

### Select Columns in Each VCF that represent a unique mutation (primary key)

### Construct a union table of distinct mutation keys (or mutation/gene keys)

### Join Genes and differential expression to mutations

### Generate Results Set






  









