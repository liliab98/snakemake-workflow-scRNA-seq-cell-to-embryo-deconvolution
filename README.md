# Snakemake workflow: single-cell RNA cell to embryo deconvolution

## Table of Contents
1. [General Info](#general-info)
2. [Usage] (#usage)
3. [Technologies](#technologies)


### General Info
***
A Snakemake workflow written by Liliane Bader as part of her Bachelor Thesis supported by Helene Kretzmer and Sandro Andreotti. The analysis steps of this workflow are part of the project [Epigenetic regulator function through mouse gastrulation]( https://github.com/HeleneKretzmer/EpigeneticRegulators_MouseGastrulation) done by [Grosswendt, S., Kretzmer, H., Smith, Z.D. et al.](https://doi.org/10.1038/s41586-020-2552-x).


## Usage

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if already available, its DOI (see above).

#### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via 

    snakemake --use-conda --cores $N

The number of threads ($N) should at least be 1, but at best 20.

#### Step 4: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push


## Technologies
***
A list of technologies used within the project:
* [STAR](https://example.com): Version 2.5.3a
* [SNPsplit](https://example.com): Version 0.3.4
* [Samtools](): Version 1.7 and 1.6
* [Perl](): Version 5.26.2
* [Bedtools](): Version 2.30.0
* [R](https://example.com): Version 4.2.0

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