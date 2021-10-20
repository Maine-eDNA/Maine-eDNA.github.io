---
layout: main
title: Metagenomics Workflow
categories: [amplicon, tutorial]
tags: [amplicon,16S,18S,metabarcoding,dada2]
permalink: /Tutorial/index
---  


{% include _amplicon_ex_workflow_toc.html %}


Here we're going to run through one way to process an amplicon dataset on a small-multi sample dataset and cover some of the many analyses for a metabarcoding pipeline.


<br>

---
---
<br>

# Pipeline Initial Assumptions
This workflow assumes that your sequencing data meets certain criteria:
- 
- 


 
<br>
<br>
<br>

---


<br>
<br>

# Working environment
If wanting to follow along, we can work on our own system if we'd like, or we can work in a "Binder" without needing to worry about setting up the appropriate environment on our own system (see next section, [Binder available](#binder-available)).

### Binder available
[Binder](https://mybinder.org/){:target="_blank"} 


<br>

## Conda setup
> **NOTE**  
> Skip this section if working in the binder environment from above. 

<br>

# Data Overview

Highland lakes info here.....

<br>
<br>


### Getting the data


Now, let's get started!
<br>
<br>

---
<br>

# Processing overview

>**Note**  
> There may be some slight differences in numbers in places due to differences in program versions between when this page was initially put together and what is in either the binder or conda environment. So don't worry if you are seeing something slightly different than what's noted or pictured anywhere below.



||Command|What it's doing|
|:--:|:--------:|----------|
|1|`cutadapt`/`filterAndTrim()`|remove primers and quality trim/filter|
|2|`learnErrors()`|generate an error model of our data|
|3|`derepFastq`|dereplicate sequences|
|4|`dada()`|infer ASVs on both forward and reverse reads independently|
|5|`mergePairs()`|merge forward and reverse reads to further refine ASVs|
|6|`makeSequenceTable()`|generate a count table|
|7|`removeBimeraDenovo()`|screen for and remove chimeras|
|8|`IdTaxa()`|assign taxonomy|


<br>
# Removing primers
To start, we need to remove the primers....


<br>
<br>
<br>
<br>
<br>


# Quality Plot Inspection
<br>
<br>

# Processing in R
<br>
### DADA2
We're going to be using [DADA2](https://benjjneb.github.io/dada2/index.html){:target="_blank"}, which is a relatively new processing workflow for recovering single-nucleotide resolved Amplicon Sequence Variants (ASVs) from amplicon data 

DADA2 R package can be found here [here](https://benjjneb.github.io/dada2/index.html){:target="_blank"}. The DADA2 team has a great tutorial available [here](https://benjjneb.github.io/dada2/tutorial.html){:target="_blank"}
<br>
<br>

<br>



### Filter and Trimming
<br>
<br>


### Dereplication
<br>
<br>

### Inferring ASVs
<br>
<br>


### Merging paired reads
<br>
<br>

### Count Table and Summary
<br>
<br>
<br>
<br>
<br>

# Assigning taxonomy
<br>
<br>
<br>
<br>
# DADA2 to Phyloseq

<br>
<br>
<br>
<br>


# Analysis in R
This portion also assumes you already have some baseline experience with R, if you aren't familiar with R at all yet it's probably a good idea to run through the [R basics page](/R/basics){:target="_blank"}. 

### Loading libraries
<br>
<br>


### Reading-In Data
<br>
<br>
<br>
<br>

# Diversity Analysis
### Beta diversity
<br>
<br>

### Alpha diversity
<br>
<br>

### Differential Abundance Analysis
<br>
<br>
<br>
<br>