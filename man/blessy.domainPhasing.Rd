% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blessy-domainPhasing.R
\name{blessy.domainPhasing}
\alias{blessy.domainPhasing}
\title{Phase Domains for Transcript-Domain Mapping Data Frame}
\usage{
blessy.domainPhasing(tx_domain_df, coordinates = TRUE)
}
\arguments{
\item{tx_domain_df}{A data frame containing transcript-domain data with the following required columns:
- \code{Transcript}: Transcript identifier.
- \code{Domain}: Domain identifier (may include NA values).
- \code{chrom}: Chromosome identifier for the domain.
- \code{chromStart}, \code{chromEnd}: Start and end positions for the domain.
- \code{strand}: Strand information for the domain.
- \code{Gene}: Gene identifier.}

\item{coordinates}{Logical flag indicating whether to include genomic coordinates and strand information
in the \code{DoCo} string. Defaults to \code{TRUE}.}
}
\value{
The input data frame with an additional column, \code{DoCo}, which contains the domain-combination
 string for each transcript.
}
\description{
This function processes a transcript-domain data frame, phasing domains along transcripts
by sorting based on strand and location to create a concatenated domain-combination (`DoCo`) string.
}
\examples{
# Example input data frame
mapping_df <- data.frame(
  Transcript = c("tx1", "tx1", "tx2"),
  Domain = c("domainA", "domainB", NA),
  chrom = c("chr1", "chr1", "chr2"),
  chromStart = c(1000, 2000, 3000),
  chromEnd = c(1100, 2100, 3100),
  strand = c("+", "-", "+"),
  Gene = c("geneA", "geneA", "geneB")
)

# Apply domain phasing with coordinates
phased_df_coords <- blessy.domainPhasing(mapping_df)

# Apply domain phasing without coordinates
phased_df_no_coords <- blessy.domainPhasing(mapping_df, coordinates = FALSE)

}
