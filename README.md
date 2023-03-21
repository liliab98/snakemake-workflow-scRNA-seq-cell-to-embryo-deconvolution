# Snakemake workflow: single-cell RNA cell to embryo deconvolution

A Snakemake workflow written by Liliane Bader as part of her Bachelor Thesis supported by Helene Kretzmer and Sandro Andreotti. The analysis steps of this workflow are part of the project ["Epigenetic regulator function through mouse gastrulation"]( https://github.com/HeleneKretzmer/EpigeneticRegulators_MouseGastrulation) following [Grosswendt, S., Kretzmer, H., Smith, Z.D. et al.](https://doi.org/10.1038/s41586-020-2552-x).

## Table of Contents
1. [Technologies](#technologies)
2. [Requirements](#requirements)
3. [Usage](#usage)


## Technologies
***
A list of technologies used within the project:
* [STAR](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html): Version 2.5.3a
* [SNPsplit](https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/SNPsplit_User_Guide.pdf): Version 0.3.4
* [Samtools](https://www.htslib.org/): Version 1.7 and 1.6
* [Perl](https://www.perl.org/): Version 5.26.2
* [Bedtools](https://bedtools.readthedocs.io/en/latest/index.html): Version 2.30.0
* [R](https://www.r-project.org/): Version 4.2.0


## Requirements
***
Multiple Samples are possible. In the following to simplify there are only two. For every sample there are two reads. They come from the same fragment and are distinguished by the attachment of  R1 or R2 to the name.


## Usage

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.

### Step 1: Obtain a copy of this workflow

If you simply want to use this workflow, download and extract the latest release. If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.
The following entries are expected in the config file:

    ID: "WT_E85"

Barcodes should be named '[ID]_barcodes.tsv.gz' and stored in the respective sample folder

    Barcodes:
        [sample1]: "resources/barcodes/[sample1]/"
        [sample2]: "resources/barcodes/[sample2]/"

Samples should be named '[sample]_R1_001.fastq.gz' and '[sample]_R2_001.fastq.gz' respectively
    
    Samples:
        [sample1]: "resources/data/[sample1]"
        [sample2]: "resources/data/[sample2]"

genome_dir is the directory where the previously generated genome indices are stored (needed for alignment)
    
    genome_dir_path: "/project/bioinf_meissner/src/SNPsplit/SNPsplit_v0.3.2" 

the following three files are needed while preparing the alignment for deconvolution
    
    all_snps_path: "resources/snps/all_SNPs_CAST_EiJ_GRCm38.txt.gz" 
    all_snps_bed_path: "resources/snps/all_SNPs_CAST_EiJ_GRCm38.bed"
    SNPwhite_list_path: "resources/scRNA/WT_SNP_white.list.tsv"

matrix_dir is where three matrices needed by the RScript "sex2embryo" for the determination of the sex
    
    matrix_dir_path: "resources/barcodes/matrix/" 

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via 

    snakemake --use-conda --cores $N

The number of threads ($N) should at least be 1, at best 20.

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/) for further details.


<!---
## Collaboration
***
Give instructions on how to collaborate with your project.
> Maybe you want to write a quote in this part. 
> It should go over several rows?
> This is how you do it.

## FAQs
***
A list of frequently asked questions
1. **This is a question in bold**
Answer of the first question with _italic words_. 
2. __Second question in bold__ 
To answer this question we use an unordered list:
* First point
* Second Point
* Third point
3. **Third question in bold**
Answer of the third question with *italic words*.
4. **Fourth question in bold**

| Headline 1 in the tablehead | Headline 2 in the tablehead | Headline 3 in the tablehead |
|:--------------|:-------------:|--------------:|
| text-align left | text-align center | text-align right |
-->