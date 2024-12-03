# longhaul

**longhaul** is a multi-purpose R package developed by the Davidson Lab. It is designed to facilitate various bioinformatics analyses focusing on long-read RNA-Seq via its component modules. 

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Overview](#overview)
- [Modules](#modules)
  - [blessy](#blessy)
    - [Introduction](#introduction)
    - [General Usage](#general-usage)
    - [Custom Usage](#custom-usage)
    - [Functions and Arguments](#functions-and-arguments)
- [License](#license)
- [Contact](#contact)

## Dependencies

The following packages are required for longhaul's modules:
- GenomicRanges
- dplyr
- tidyr
- UCSC.utils
- stats

Which can be installed via:

```R
# Install core dependencies from CRAN
install.packages(c("dplyr", "tidyr"))

# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "UCSC.utils"))

```

## Installation

You can install the **longhaul** package from GitHub using the `devtools` package:

```R
# Install devtools if you haven't already
install.packages("devtools")

# Install longhaul from GitHub
devtools::install_github("DavidsonGroup/longhaul")
```

## Overview

**longhaul** is a comprehensive toolkit intended to streamline various computational biology workflows within the Davidson Lab. It encompasses multiple modules, each tailored to specific tasks in long-read RNA-Seq analysis.

The **longhaul** package includes the following modules:

- blessy: a tool for differential analysis on phased protein domains.


## Modules

### blessy 

#### Introduction
blessy is a specialized module within the longhaul package designed for performing differential analysis on a novel genomic feature termed 'Domain Combination'
or 'DoCo' for short. 

*pic of Domain Combination*


#### General use

#### Custom use



#### Functions and Arguments:
- **Fetch Annotation Tracks**: Retrieve transcript and domain annotation tracks from the UCSC Genome Browser.
- **Custom Annotations**: Use custom BED-like data frames for transcript and domain annotations.
- **Domain Mapping**: Map domains to transcripts based on genomic overlaps.
- **Deduplication**: Remove non-exact domain mappings to ensure accuracy.
- **Phasing**: Phase domains along transcripts to create Domain Combination (DoCo) strings.
- **Aggregation**: Generate DoCo-level counts from transcript-level RNA-seq counts.
- **Phasing Dictionary**: Create a comprehensive dictionary summarizing gene, DoCo, and transcript relationships.



