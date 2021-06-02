# Expression and SNPs to clinic
Smooth transition of called variants from RNAseq/DNAseq and expression to the clinic. 

## Contributors 

- Jyoti Kataria - Sysadmin 
- Amit Yadav - Sysadmin 
- Adam Pater Faranda - Data Guru, Tech support 
- Katarina - Data Guru 
- Varuna Chander - Writer, Tech support 
- Kym Pagel - Lead, Liaison 

## Goal 
Build a streamlined and easy to use workflow for reporting expressed variants from RNAseq to the clinic. 


## Introduction 


We use the CTAT-Mutation pipeline will be used to call expresseed variants from RNAseq data. The CTAT-Mutation pipeline (https://github.com/NCIP/ctat-mutations/wiki) makes it easy to discover variants from RNA-seq data, and requires only the RNA-seq reads as input. The pipeline also annotates variants, including the RADAR and RediPortal databases for identifying likely RNA-editing events, dbSNP and gnomAD for annotating common variants, COSMIC to highlight known cancer mutations, and OpenCRAVAT to annotate and prioritize variants according to likely biological impact and relevance to cancer. The CTAT-Mutations Pipeline integrates GATK Best Practices along with downstream steps to annotate and filter variants, and to additionally prioritize variants that may be relevant to cancer biology. 

We will then use the GATK Best Practices pipeline to call variants from DNAseq. Next, we will identify genes that are differentially expressed. Finally, we will aggregate the variants identified through DNA and RNAseq, and curate extensive clinical annotations using OpenCRAVAT to identify priority variants. 

## Test Data 

In order to thoroughly evaluate this the pipeline we've identified several studies that meet the following basic criteria
- the study provides a matching DNA Sequencing and RNA Seqeuncing reads from same biological sample
- the study provides a meaningful contrast for evaluating differential expression
- reads are paired end.

We've identified a study archived in the Gene Expression Omnibus (Geo Accession:GSE75935) that provides RNA Seq and DNA Seq data from human tissue samples collected from ovarian cancer patients. The data includess three cancer patients, three tumor samples per patient from different sites, two normal tissue samples from two different patients, four cell lines.

## Installation 


## Methods
1. Obtain high quality test data
1. Construct the variant calling and annotation pipeline
1. Construct the gene expression analysis pipeline 
1. Document the pipeline thoroughly 
1. Construct output file for OMOP group

### Implementation 



#### Inputs 

DNAseq Workflow:

BWA FASTA indexer: Indexes a FASTA reference genome sequence for downstream use by the BWA-backtrack and BWA-MEM mapping apps.
Input: Reference FASTA file (UCSC hg19)
Output:BWA FASTA index file 

BWA-MEM mapping: Maps FASTQ (paired or unpaired reads) to reference genome using BWA-MEM algorithm.

Inputs:
- DNA paired end fastq files
- BWA FASTA index file 

Outputs: 
- Sorted BAM file
- Index BAM file


Exome GATK lite pipeline: Deduplicates, realigns and recalibrates base quality scores, and calls SNPs and indels in the human exome. Runs the following software from the Picard 1.104 and GATK-lite v2.3 suite of tools:
- Picard MarkDuplicates
- GATK-lite RealignerTargetCreator (if configured to do realignment around novel indels)
- GATK-lite IndelRealigner
- GATK-lite BaseRecalibrator
- GATK-lite PrintReads
- GATK-lite UnifiedGenotyper

Input:
- Requires a sorted BAM file.

Outputs:
- gzipped VCF file with the called variants
- Variants indexed file

----------//-------------//---------------------

RNAseq Workflow:

I. Trinity CTAT:
- RNAseq fastq file
- Reference genome


OpenCRAVAT



II. DESeq2:
- RNAseq fastq file


ViraVate2



CombineR



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
<img width="303" alt="Screen Shot 2021-06-02 at 6 24 51 PM" src="https://user-images.githubusercontent.com/5508556/120559925-e4d8eb80-c3cf-11eb-9d36-835f990482f2.png">
Flowchart of the pipeline

## Results 


## References 

- Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75935 
- CTAT pipeline https://github.com/NCIP/ctat-mutations/wiki 
- OpenCRAVAT https://www.cancergeneticsjournal.org/article/S2210-7762(20)30193-9/pdf 
- GATK Best Practices https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows 
- DSeq2 https://bioconductor.org/packages/release/bioc/html/DESeq2.html 
- ViraVate https://github.com/collaborativebioinformatics/Differential_Expression_and_Variant_Association 
- DNANexus documentation https://documentation.dnanexus.com/developer/apps/execution-environment/connecting-to-jobs 
