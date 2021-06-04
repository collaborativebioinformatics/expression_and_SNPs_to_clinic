# File: Combine.R
# Purpose: Integrate results from Variant Calling and Differential Expression
#          workflows into a unified tsv file for each sample.
#          
# Usage:   Update path, set infile names to catch pipeline outputs
# Author:  Adam Faranda
# Created: June 3, 2021
###################################################################################

## Update Github:

## Filter For PASS, No QUAL
## One DNA, 2xRNA use RNA_1 & RNA_2
## USE Concatenated Chromosome-position-ref-alt as unique mutation ID
## Use Gene SYmbol to map Diff Exp, greatest sig fold change to dedup
## Use the columns for each type of HGVS c, p
## CRV fields are fixed position
## Transcript Field (ENST)-- from MANE Canonical
## COSMIC, Clinvar, 


## Grab top 5 DE Genes -- return separate output file
## Gene Name, Fold Change, unique mutation count 






###################### Setup Environment and Import Tables ########################
if(!"dplyr" %in% row.names(installed.packages())){
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if(!"tidyr" %in% row.names(installed.packages())){
  install.packages("tidyr", repos = "https://cloud.r-project.org")
}
if(!"tibble" %in% row.names(installed.packages())){
  install.packages("tibble", repos = "https://cloud.r-project.org")
}

library(dplyr)
library(plyranges)

wd <- "~/expression_and_SNPs_to_clinic/"   # Update Path depending on system
setwd(wd)

### GATK Variants from DNA Seq
dna_fn <- "OCinput.vcf.gz.vcf"
system(
  paste(wd,"clean_vcf.sh ", dna_fn, sep="")
)
dna_seq <- read.table(
  paste0(wd,dna_fn,"_Clean.tsv"),
  quote="", header=F, sep="\t"
)

### GATK Variants from RNA Seq File 1
rna_fn_1 <- "rna_1_vcf.vcf"
system(
  paste(wd,"clean_vcf.sh ", rna_fn_1, sep="")
)
rna_seq_1 <- read.table(
  paste0(wd,rna_fn_1,"_Clean.tsv"),
  quote="", header=F, sep="\t"
)

### GATK Variants from RNA Seq File 1
rna_fn_2 <- "rna_2_vcf.vcf"
system(
  paste(wd,"clean_vcf.sh ", rna_fn_2, sep="")
)
rna_seq_2 <- read.table(
  paste0(wd,rna_fn_2,"_Clean.tsv"),
  quote="", header=F, sep="\t"
)


### Differential Expression (Tumor vs Normal)
deg_fn <- "Group1_vs_Group2.allGenes.csv"
deg_table <- read.table(deg_fn, sep="\t", header=T, quote="")

######################### remove after testing ##################
rna_seq_1 <- dna_seq[1:15,]
rna_seq_2 <- dna_seq[5:20,]
dna_seq <- dna_seq[10:25,]
########################### Cleanup Column Headers ###########################
var_table_headers <- function(var_table){
  var_table %>%
    select(
      CHROM=V6,
      POS=V10,
      REF=V11,
      ALT=V12,
      HUGO_SYMBOL=V1,
      TRANSCRIPT=V2,
      SEQ_ONT=V3,
      TRANSCRIPT_HGVS=V4,
      PROTEIN_HVGS=V5
    ) %>%
    mutate(
      VAR_ID=paste0(CHROM,":",POS,"_",REF,"/",ALT)
    )
}

dna_seq_check <- var_table_headers(dna_seq) %>%
  filter(
    grepl("^[0-9].*$", POS)
  ) 
rna_seq_1 <- var_table_headers(rna_seq_1)
rna_seq_2 <- var_table_headers(rna_seq_2)
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

## CRV String Elements


### Select Columns in Each VCF that represent a unique mutation (primary key)
### Construct a union table of distinct mutation keys (or mutation/gene keys)
## Collapse to distinct variants and use "Source" to
## indicate which pipeline detected the variation (DNA, RNA, BOTH)
variants <- bind_rows(
  dna_seq %>%
    mutate(Source = "DNA"),
  rna_seq_1 %>%
    mutate(Source = "RNA_1"),
  rna_seq_2 %>%
    mutate(Source = "RNA_2")
) %>%
  group_by(
    mutation_id, REF, ALT, QUAL, FILTER
  ) %>%
  summarise(
    Source=ifelse(n()==1, Source, "BOTH")
  ) %>% as.data.frame()


### Join Genes and differential expression to mutations

### Generate Results Set


################ MOCKUP FOR SCREENSHOTS ###########################
dna_seq_check %>%
  filter(HUGO_SYMBOL=="SAMD11") %>%
  mutate(
    OBSERVED_IN_DNA=c(
      "YES", "YES", "NO", "NO", "NO", "YES", "NO", "NO"
    ),
    OBSERVED_IN_RNA_1=c(
      "YES", "NO", "NO", "YES", "NO", "YES", "NO", "NO"
    ),
    OBSERVED_IN_RNA_2=c(
      "YES", "NO", "NO", "NO", "NO", "YES", "NO", "NO"
    ),
    COSMIC_ID=c("COSV52316050"),
    ZYGOSITY=c("het", "het", "hom", "hom", "hom", "hom", "hom", "hom"),
    DE_LOG_FC=10.5,
    DE_FDR=0.05
  )

top_five=data.frame(
  HUGO_SYMBOL=c("RUNX1", "SAMD11", "TP53", "ATF4", "EGR1"),
  FOLD_CHANGE=c(12, 10.5, -8, 6, 5),
  MUTATION_COUNT=c(3,8, 15, 2, 1)
)




