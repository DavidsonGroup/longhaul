#' Map Domains to Transcripts By Finding GRangesLists' Overlaps
#'
#' This function maps domains to transcripts by identifying overlaps between two GRangesLists representing transcripts and domains.
#' It returns a mapping data frame containing all columns from both inputs for matching rows.
#' Column names are renamed for clarity and downstream analysis, and duplicate columns are removed.
#'
#' @param tx_grangesList A `GRangesList` or `CompressedGRangesList` object representing transcripts.
#' @param domain_grangesList A `GRangesList` or `CompressedGRangesList` object representing domains.
#' @param tx_df A BED-like data frame representing transcripts, including columns like `chrom`, `chromStart`, `chromEnd`, `strand`, etc.
#' @param domain_df A BED-like data frame representing domains, including columns like `chrom`, `chromStart`, `chromEnd`, `strand`, etc.
#'
#' @return A data frame with where each row represents mapping of a domain to a transcripts.
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' 
#' @examples
#' # Example BED-like data frames
#' tx_df <- data.frame(
#'   chrom = c("chr1", "chr1", "chr2"),
#'   chromStart = c(1000, 2000, 3000),
#'   chromEnd = c(1500, 2500, 3500),
#'   name = c("tx1", "tx2", "tx3"),
#'   strand = c("+", "-", "+"),
#'   thickStart = c(1000, 2000, 3000),
#'   thickEnd = c(1500, 2500, 3500),
#'   blockCount = c(2, 2, 2),
#'   blockSizes = c("100,200", "150,250", "200,300"),
#'   blockStarts = c("0,400", "0,500", "0,600"),
#'   geneName = c("geneA", "geneB", "geneC"),
#'   stringsAsFactors = FALSE
#' )
#'
#' domain_df <- data.frame(
#'   chrom = c("chr1", "chr1", "chr2"),
#'   chromStart = c(1200, 2100, 3100),
#'   chromEnd = c(1300, 2200, 3200),
#'   name = c("domainA", "domainB", "domainC"),
#'   strand = c("+", "-", "+"),
#'   blockStarts = c("0", "0", "0"),
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Example GRangesLists
#' tx_grangesList <- blessy.dfToGRangesList(tx_df)
#' 
#' domain_grangesList <- blessy.dfToGRangesList(domain_df)
#' 
#' # Map domains to transcripts and combine data frames
#' intersection_df <- blessy.mapDomainToTranscript(tx_grangesList, domain_grangesList, tx_df, domain_df)
#'
#' @export
blessy.mapDomainToTranscript <- function(tx_grangesList, domain_grangesList, tx_df, domain_df) {
  
  # Find overlaps between the two GRangesLists
  overlaps <- findOverlaps(tx_grangesList, domain_grangesList)
  
  # Extract indices of overlapping elements
  query_hits <- queryHits(overlaps)     # Indices in tx_grangesListlist (tx_df)
  subject_hits <- subjectHits(overlaps) # Indices in domain_grangesListlist (domain_df)
  
  # Extract matching rows from tx_df and domain_df using the indices
  matched_tx <- tx_df[query_hits, , drop = FALSE]
  matched_domain <- domain_df[subject_hits, , drop = FALSE]
  
  # Rename columns for the first data frame (transcripts)
  tx_rename <- c(
    "chromStart" = "txStart",
    "chromEnd" = "txEnd",
    "name" = "Transcript",
    "thickStart" = "cdsStart",
    "thickEnd" = "cdsEnd",
    "blockCount" = "exonCount",
    "blockSizes" = "exonSizes",
    "blockStarts" = "exonRelativeStarts",
    "geneName" = "Gene"
  )
  colnames(matched_tx) <- ifelse(colnames(matched_tx) %in% names(tx_rename), 
                                 tx_rename[colnames(matched_tx)], colnames(matched_tx))
  
  # Rename columns for the second data frame (domains)
  domain_rename <- c(
    "name" = "Domain",
    "blockStarts" = "chromStarts"
  )
  colnames(matched_domain) <- ifelse(colnames(matched_domain) %in% names(domain_rename), 
                                     domain_rename[colnames(matched_domain)], colnames(matched_domain))
  
  # Combine the matched data frames into one
  combined_df <- cbind(matched_tx, matched_domain)
  
  # Remove duplicated columns (not renamed)
  duplicated_columns <- intersect(names(matched_tx), names(matched_domain))
  combined_df <- combined_df[, !duplicated(names(combined_df))]
  
  return(combined_df)
}

