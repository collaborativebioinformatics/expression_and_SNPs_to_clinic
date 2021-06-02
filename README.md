# Expression and SNPs to clinic
Smooth transition of called variants from RNAseq/DNAseq and expression to the clinic. 

## Contributors 

- Jyoti Kataria 
- Amit Yadav
- Adam Pater Faranda
- Katarina
- Varuna Chander 
- Kym Pagel 

## Goal 
Build a streamlined and easy to use workflow for reporting expressed variants from RNAseq to the clinic. 


## Introduction 


We use the CTAT-Mutation pipeline will be used to call expresseed variants from RNAseq data. The CTAT-Mutation pipeline (https://github.com/NCIP/ctat-mutations/wiki) makes it easy to discover variants from RNA-seq data, and requires only the RNA-seq reads as input. The pipeline also annotates variants, including the RADAR and RediPortal databases for identifying likely RNA-editing events, dbSNP and gnomAD for annotating common variants, COSMIC to highlight known cancer mutations, and OpenCRAVAT to annotate and prioritize variants according to likely biological impact and relevance to cancer. The CTAT-Mutations Pipeline integrates GATK Best Practices along with downstream steps to annotate and filter variants, and to additionally prioritize variants that may be relevant to cancer biology. 

We will then use the GATK Best Practices pipeline to call variants from DNAseq. Next, we will identify genes that are differentially expressed. Finally, we will aggregate the variants identified through DNA and RNAseq, and curate extensive clinical annotations using OpenCRAVAT to identify priority variants. 

## Installation 


## Methods
1. Obtain high quality test data
1. Construct the variant calling and annotation pipeline
1. Construct the gene expression analysis pipeline 
1. Document the pipeline thoroughly 
1. Construct output file for OMOP group

### Implementation 

#### Inputs 

GATK:
- DNAseq fastq file
- Reference genome


Trinity CTAT:
- RNAseq fastq file
- Reference genome


Dseq2:
RNAseq fastq file

#### Outputs 

A TSV-delimited file per sample. Each line describes one variant, including the following fields: 
- Variant HGVS expression 
- Hugo gene identifier 
- Variant sequence ontology consequence
- ClinVar annotation
- COSMIC annotation 
- Variant zygosity 
- Gene-level differential expression value
- Source (DNA,RNA,Both)
- Source identifier / Tissue origin, if relevant 

## Operation 


## Flowchart
<img width="448" alt="flowchart" src="https://user-images.githubusercontent.com/5508556/120391842-a1ae4800-c2fd-11eb-92fe-06baed291a44.png">
Flowchart of the pipeline

## Results 


## References 

- Data ??? 
- CTAT pipeline https://github.com/NCIP/ctat-mutations/wiki 
- OpenCRAVAT https://www.cancergeneticsjournal.org/article/S2210-7762(20)30193-9/pdf 
- GATK Best Practices https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows 
- DSeq2 https://bioconductor.org/packages/release/bioc/html/DESeq2.html 
- ViraVate https://github.com/collaborativebioinformatics/Differential_Expression_and_Variant_Association 
- DNANexus documentation https://documentation.dnanexus.com/developer/apps/execution-environment/connecting-to-jobs 
