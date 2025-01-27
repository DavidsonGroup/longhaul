#' Map Domains to Transcripts Using BED-Like Data Frames
#'
#' This function takes two BED-like data frames (one for transcripts, one for domains),
#' converts them internally to `GRangesList` objects, then finds overlaps to map domains
#' to transcripts. It returns a merged data frame with renamed columns for clarity.
#'
#' @param tx_df A BED-like data frame representing transcripts. **Required** columns:
#'   - \code{chrom}: Chromosome identifier.
#'   - \code{thickStart}, \code{thickEnd}: Genomic start and end coordinates used for overlap.
#'   - \code{strand}: Strand information (+, -, or *).
#'   - \code{exonStarts}, \code{exonEnds} (for downstream domain deduplication).
#'   - \code{blockStarts}, \code{blockEnds} (for downstream domain deduplication).
#'   Optionally, additional fields (e.g., \code{chromStart}, \code{chromEnd}, \code{blockCount}, 
#'   \code{name}, \code{geneName}) may be included and will be renamed or preserved as needed.
#'
#' @param domain_df A BED-like data frame representing domains. **Required** columns:
#'   - \code{chrom}: Chromosome identifier.
#'   - \code{thickStart}, \code{thickEnd}: Genomic start and end coordinates used for overlap.
#'   - \code{strand}: Strand information (+, -, or *).
#'   - \code{blockStarts}, \code{blockEnds}: Comma-separated domain block start/end coordinates
#'     (needed for downstream domain deduplication).
#'   Optionally, additional fields (e.g., \code{name}) may be included.
#'
#' @return A data frame in which each row represents a mapped domain-to-transcript overlap.
#'         Columns from both \code{tx_df} and \code{domain_df} are combined; some are renamed
#'         for clarity and downstream processing (e.g., domain deduplication).
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @export
blessy.mapDomainToTranscript <- function(tx_df, domain_df) {
  
  # 1) Convert data frames to GRangesList objects using your existing helper
  tx_grangesList <- blessy.dfToGRangesList(tx_df)
  domain_grangesList <- blessy.dfToGRangesList(domain_df)
  
  # 2) Find overlaps between the GRangesList objects
  overlaps <- findOverlaps(tx_grangesList, domain_grangesList)
  
  # 3) Extract indices of overlapping elements
  query_hits <- queryHits(overlaps)     # Indices in tx_grangesList (matching rows in tx_df)
  subject_hits <- subjectHits(overlaps) # Indices in domain_grangesList (matching rows in domain_df)
  
  # 4) Match rows from tx_df and domain_df
  matched_tx <- tx_df[query_hits, , drop = FALSE]
  matched_domain <- domain_df[subject_hits, , drop = FALSE]
  
  # 5) Optionally rename columns for clarity
  tx_rename <- c(
    "chromStart"   = "txStart",
    "chromEnd"     = "txEnd",
    "name"         = "Transcript",
    "thickStart"   = "cdsStart",
    "thickEnd"     = "cdsEnd",
    "blockCount"   = "exonCount",
    "blockSizes"   = "exonSizes",
    "blockStarts"  = "exonRelativeStarts",
    "geneName"     = "Gene"
    # etc.
  )
  colnames(matched_tx) <- ifelse(
    colnames(matched_tx) %in% names(tx_rename),
    tx_rename[colnames(matched_tx)],
    colnames(matched_tx)
  )
  
  domain_rename <- c(
    "name"         = "Domain",
    "blockStarts"  = "domainBlockStarts",
    "blockEnds"    = "domainBlockEnds"
    # etc.
  )
  colnames(matched_domain) <- ifelse(
    colnames(matched_domain) %in% names(domain_rename),
    domain_rename[colnames(matched_domain)],
    colnames(matched_domain)
  )
  
  # 6) Combine the matched data frames
  combined_df <- cbind(matched_tx, matched_domain)
  
  # 7) Remove duplicate columns
  duplicated_columns <- intersect(names(matched_tx), names(matched_domain))
  combined_df <- combined_df[, !duplicated(names(combined_df))]
  
  return(combined_df)
}
