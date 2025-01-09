#' Domain Deduplication for Transcript-Domain Mapping
#'
#' This function deduplicates a transcript-domain mapping data frame by filtering out 
#' domains with non-exact matching. It ensures that only matches with relevant
#' exon-block relationships are retained.
#'
#' @param tx_domain_df A data frame containing transcript-domain mapping data, including columns:
#'   - \code{exonStarts}, \code{exonEnds}: Comma-separated exon start and end positions.
#'   - \code{blockStarts}, \code{blockEnds}: Comma-separated block start and end positions.
#' @param unique_domain A logical value indicating whether to further deduplicate domains by keeping only unique domains per transcript. Default is \code{FALSE}.
#'
#' @return A deduplicated data frame of transcript-domain mapping.
#'
#' @examples
#' # Example data frame
#' intersect_df <- data.frame(
#'   exonStarts = c("1000,2000,3000", "1500,2500,3500"),
#'   exonEnds = c("1100,2100,3100", "1600,2600,3600"),
#'   blockStarts = c("1000,2000,3000", "1500,2500,3500"),
#'   blockEnds = c("1100,2100,3100", "1600,2600,3600")
#' )
#'
#' # Apply the function to deduplicate domains
#' refined_df <- blessy.domainDeduplication(intersect_df)
#'
#' @export
blessy.domainDeduplication <- function(tx_domain_df, unique_domain = FALSE) {
  # Parse coordinates into lists of numeric vectors
  exon_starts_list <- strsplit(as.character(tx_domain_df$exonStarts), ",")
  exon_ends_list <- strsplit(as.character(tx_domain_df$exonEnds), ",")
  block_starts_list <- strsplit(as.character(tx_domain_df$blockStarts), ",")
  block_ends_list <- strsplit(as.character(tx_domain_df$blockEnds), ",")
  
  # Convert all to numeric at once
  exon_starts_list <- lapply(exon_starts_list, as.numeric)
  exon_ends_list <- lapply(exon_ends_list, as.numeric)
  block_starts_list <- lapply(block_starts_list, as.numeric)
  block_ends_list <- lapply(block_ends_list, as.numeric)
  
  is_valid_domain <- mapply(function(exon_starts, exon_ends, block_starts, block_ends) {
    # Number of blocks and exons
    n_blocks <- length(block_starts)
    n_exons <- length(exon_starts)
    
    # Helper function to find which exons a block could belong to
    # A block belongs to an exon if blockStart >= exonStart and blockEnd <= exonEnd
    find_exons_for_block <- function(bstart, bend) {
      which(bstart >= exon_starts & bend <= exon_ends)
    }
    
    # Single-block domain validation
    if (n_blocks == 1) {
      candidates <- find_exons_for_block(block_starts, block_ends)
      # Must map exactly to one exon
      return(length(candidates) == 1)
    }
    
    # Multi-block domain validation:
    # Each block should map exactly to one exon with specific conditions:
    # - First block: blockEnds == exonEnds
    # - Last block: blockStarts == exonStarts
    # - Intermediate blocks: exact match of both start and end
    
    mapped_exons <- integer(n_blocks)  # store the exon index each block maps to
    
    for (i in seq_len(n_blocks)) {
      candidates <- find_exons_for_block(block_starts[i], block_ends[i])
      
      # If no or multiple candidate exons, invalid
      if (length(candidates) != 1) {
        return(FALSE)
      }
      
      # Check boundary conditions
      exon_idx <- candidates
      if (i == 1) {
        # First block: must end at exon boundary
        if (block_ends[i] != exon_ends[exon_idx]) return(FALSE)
      } else if (i == n_blocks) {
        # Last block: must start at exon boundary
        if (block_starts[i] != exon_starts[exon_idx]) return(FALSE)
      } else {
        # Intermediate blocks: must match exon start/end exactly
        if (block_starts[i] != exon_starts[exon_idx] || block_ends[i] != exon_ends[exon_idx]) return(FALSE)
      }
      
      mapped_exons[i] <- exon_idx
    }
    
    # Check that exons mapped by consecutive blocks are consecutive integers (no exon skipping)
    # This means the differences between consecutive mapped exon indices should be exactly 1
    if (!all(diff(mapped_exons) == 1)) return(FALSE)
    
    # If all conditions are met
    return(TRUE)
  }, exon_starts_list, exon_ends_list, block_starts_list, block_ends_list)
  
  # Filter the data frame to keep only valid domains
  tx_domain_df <- tx_domain_df[is_valid_domain, , drop = FALSE]
  
  # Deduplicate rows by Transcript and Domain if unique_domain is TRUE
  if (unique_domain) {
    tx_domain_df <- tx_domain_df[!duplicated(tx_domain_df[, c("Transcript", "Domain")]), ]
  }
  
  return(tx_domain_df)
}
