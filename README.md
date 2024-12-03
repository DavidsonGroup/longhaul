# longhaul

**longhaul** is a multi-purpose R package developed by the Davidson Lab. It is designed to facilitate various bioinformatics analyses focusing on long-read RNA-Seq via its component modules. 

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Overview](#overview)
- [Modules](#modules)
  - [blessy](#blessy)
    - [General Usage](#general-usage)
        - [Using pre-defined UCSC annotation tracks](#using-pre-defined-ucsc-annotation-tracks)
        - [Using custom annotation tracks](#using-custom-annotation-tracks)
        - [Output](#output)
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

*blessy* is a module within the longhaul package for performing differential analysis on a novel genomic feature termed 'Domain Combination'
or 'DoCo' for short. 

In the figure below, the hypothetical gene A has 4 transcripts (Tx1, Tx2, Tx3, Tx4). On each transcript, regions that encode for protein domains (PD1, PD2) can be found. Transcripts with the similar domain phasing are categorized into the same 'Domain Combination' or 'DoCo' group. DoCo class is then considered per gene, making DoCo an intermediate feature of gene and transcript. Of note, transcripts with no domain (Tx4) are still categorized into an ‘empty’ DoCo group (;;; geneA).

The DoCo class can be used to group transcripts of an existing RNA-Seq count. Besides making biologically equivalent transcript groups, the count of each DoCo is the aggregated count of its component transcripts, addressing the previous issue of low count per feature in differential transcript analyses such as DTE or DTU. 

![The Concept of Domain Combination](figures/DoCo_concept.png)


#### General Usage
*blessy* requires annotation tracks for transcripts and protein domains to create the DoCo class. Additionally, an RNA-Seq transcript count must be provided for aggregating count at transcript-level to DoCo-level, once the information on DoCo class is generated. 

##### Using pre-defined UCSC annotation tracks
The default mode to run *blessy* is using pre-defined annotation tracks on the UCSC Genome Browser and an user-provided transcript count. *blessy* will return a list containing two R data frames: a dictionary showing the hierarchical relationship of gene, DoCo and transcript from the annotations of choice, and the count table at DoCo level gained from aggregating the transcript count using this dictionary. If you would like to know which annotations can be used for blessy, you can find them [HERE](https://genome.ucsc.edu/cgi-bin/hgTables). We highly recommend using the GENCODE or NCBI RefSeq tracks for transcript annotation and UniProt or Pfam tracks for domain annotation, as *blessy* is tailored around these annotations.


```R
blessy_outs <- blessy(genomeAssembly, transcriptAnnotation, domainAnnotation, transcriptCount)
```

**genomeAssembly** - A string specifying the UCSC assembly identifier used for the two transcriptAnnotation and domainAnnotation below. This string corresponds to the identifier which can be found at the end of the 'Assembly' option in the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables), or the 'Database' seen in the 'Data format description' page of a specific track. (e.g., "hg38").

**transcriptAnnotation** - A string specifying the UCSC table identifier of a transcript track. This string corresponds to the 'Table' option in the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables), or the 'Primary Table' seen in the 'Data format description' page of a specific track. (e.g., "wgEncodeGencodeBasicV44").

**domainAnnotation** - A string specifying the UCSC table identifier of a domain track. This string corresponds to the 'Table' option in the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables), or the 'Primary Table' seen in the 'Data format description' page of a specific track. (e.g., "unipDomain").

**transcriptCount** - An R data frame containing RNA-Seq transcript count. The first column must be named 'TranscriptID' storing string values of transcript identifiers. Other columns are considered numeric count across different biological samples.

##### Using custom annotation tracks
Here *blessy* allows users to provide their own annotation tracks of choice, leaving room for annotation customization. For example, users can choose to retrieve multiple domain annotation tracks via blessy.getDomainTrack() (see [Functions and Arguments](#functions-and-arguments)), and integrate these domain tracks into one comprehensive domain annotation using simple R commands:

```R
# Fetch domain tracks on the UCSC Genome Browser
unipDomain <- blessy.getDomainTrack("hg38", "unipDomain")
unipInterest <- blessy.getDomainTrack("hg38", "unipInterest")

# Combine the two data frames
combined_df <- bind_rows(unipDomain, unipInterest)

# Filter out rows containing "Disordered"
filtered_df <- combined_df %>% 
  filter(!grepl("Disordered", name)) 

# Sort 
domain_df <- filtered_df %>% 
  arrange(chrom, chromStart) 
```

To run *blessy* using annotations not from UCSC, users must ensure that these annotations are first read into BED-like R data frames. These data frames must include the following columns:

**chrom**- Chromosome or scaffold name of the feature.

**chromStart** - Start coordinate on the chromosome or scaffold for the feature considered.

**chromEnd** - End coordinate on the chromosome or scaffold for the feature considered.

**blockCount** - Number of blocks (exon for transcripts, domain-block for domains) on the line of the BED-like data frame.

**blockSizes** - List of values separated by commas corresponding to the size of the blocks (the number of values must correspond to that of the 'blockCount').

**blockStarts** - List of values separated by commas corresponding to the starting coordinates of the blocks, coordinates calculated relative to those present in the chromStart column (the number of values must correspond to that of the 'blockCount').

**thickStart** - **(for transcript only)** Start coordinate of the coding sequence of the transcript.

**thickEnd** - **(for transcript only)** End coordinate of the coding sequence of the transcript.

**geneName** - **(for transcript only)** Identifier of the gene the transcript belongs to.


*blessy* can then use these annotations to define the DoCo class and aggregate the provided transcript count:

```R
blessy_outs_custom <- blessy.usingCustomAnnotation(customTranscriptAnnotation, customDomainAnnotation, transcriptCount)
```

**customTranscriptAnnotation** - A BED-like data frame representing transcript annotation of the required format. 

**customDomainAnnotation** - A BED-like data frame representing domain annotation of the required format. 

**transcriptCount** - An R data frame containing RNA-Seq transcript count. The first column must be named 'TranscriptID' storing string values of transcript identifiers. Other columns are considered numeric count across different biological samples.


##### Output

*blessy* returns a list containing two data frames: 
  - **phasing_dict** - A dictionary showing the hierarchical relationship of gene, DoCo and transcript from the annotations of choice
  - **doco_count** - A count table at DoCo level

#### Functions and Arguments:
Below we outline the component functions of blessy along with use-cases for each, to assist users in customizing the module to fit their specific needs.

- **Fetch Annotation Tracks**: Retrieve transcript and domain annotation tracks from the UCSC Genome Browser.
- **Custom Annotations**: Use custom BED-like data frames for transcript and domain annotations.
- **Domain Mapping**: Map domains to transcripts based on genomic overlaps.
- **Deduplication**: Remove non-exact domain mappings to ensure accuracy.
- **Phasing**: Phase domains along transcripts to create Domain Combination (DoCo) strings.
- **Aggregation**: Generate DoCo-level counts from transcript-level RNA-seq counts.
- **Phasing Dictionary**: Create a comprehensive dictionary summarizing gene, DoCo, and transcript relationships.



