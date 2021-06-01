# Expression and SNPs to clinic
Smooth transition of called variants from RNAseq and expression to the clinic

## Contributors 

- Jyoti Kataria 
- Amit Yadav
- Adam Pater Faranda
- Katarina
- Varuna Chander 
- Kym Pagel 

## Goal 
Build a streamlined and easy to use workflow for reporting expressed variants from RNAseq to the clinic


## Introduction 


## Installation 


## Methods
- Choice of Test data
- Variant Calling and Annotation Pipeline
- Gene expression analysis
- Identification of variants with clinical relevance
- Visualization of data

### Implementation

#### Inputs 


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

