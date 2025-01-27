#' blessy's Wrapper Function using Custom Annotation
#'
#' This wrapper function performs the workflow starting from mapping domains to transcripts
#' using user-provided transcript and domain annotation data frames of BED-like format or GRangesList objects.
#' It creates a phasing dictionary and generates DoCo-level count data frames.
#'
#' @param tx_df A BED-like data frame representing transcript annotation or a GRangesList object. The annotation must include information on exon locations in \code{blockCount}, \code{blockSizes}, \code{blockStarts}, and gene name in \code{geneName}.
#' @param domain_df A BED-like data frame representing domain annotation or a GRangesList object. The annotation must include information on domain block locations in \code{blockCount}, \code{blockSizes}, \code{blockStarts}.
#' @param tx_count A data frame containing transcript counts with:
#'   - \code{TranscriptID}: Transcript IDs (required as the first column).
#'   - Additional columns representing RNA-seq counts for biological samples (must be numeric).
#' @param unique_domain Logical flag indicating whether to deduplicate domains by keeping only unique domains per transcript. 
#'   Defaults to \code{FALSE}.
#' @param coordinates Logical flag indicating whether to include genomic coordinates and strand information
#'   in the \code{DoCo} string during the phasing step. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#'   - \code{phasing_dict}: The phasing dictionary data frame.
#'   - \code{doco_count}: The DoCo-level count data frame.
#'
#' @export
blessy.usingCustomAnnotation <- function(tx_df, domain_df, tx_count, unique_domain = FALSE, coordinates = TRUE) {
  # Validation for tx_df
  cat("Step 1/8: Validating transcript data frame format...\n")
  tx_required_cols <- c("chrom", "txStart", "txEnd", "strand", "blockCount", "blockSizes", "blockStarts", "geneName")
  if (!all(tx_required_cols %in% colnames(tx_df))) {
    stop("The transcript data frame (tx_df) must contain the following columns: ", paste(tx_required_cols, collapse = ", "))
  }
  
  # Validation for domain_df
  cat("Step 2/8: Validating domain data frame format...\n")
  domain_required_cols <- c("chrom", "chromStart", "chromEnd", "strand", "blockCount", "blockSizes", "blockStarts")
  if (!all(domain_required_cols %in% colnames(domain_df))) {
    stop("The domain data frame (domain_df) must contain the following columns: ", paste(domain_required_cols, collapse = ", "))
  }
  
  # Step 3: Map domains to transcripts
  cat("Step 3/8: Matching domains to transcripts...\n")
  mapped_df <- blessy.mapDomainToTranscript(tx_df, domain_df)
  
  # Step 4: Add exon and block starts/ends
  cat("Step 4/8: Adding exon and block starts/ends...\n")
  cat("Note: Patience is bitter, but its fruit is sweet. \n")
  starts_ends_df <- blessy.addStartsEnds(mapped_df)
  
  # Step 5: Deduplicate domain mappings
  cat("Step 5/8: Deduplicating domain mappings...\n")
  deduplicated_df <- blessy.domainDeduplication(starts_ends_df, unique_domain = unique_domain)
  
  # Step 6: Create phasing information with the 'coordinates' parameter
  cat("Step 6/8: Creating phasing information...\n", coordinates)
  phased_df <- blessy.domainPhasing(deduplicated_df, coordinates = coordinates)
  
  # Step 7: Create the phasing dictionary
  cat("Step 7/8: Creating the phasing dictionary...\n")
  phasing_dict <- blessy.createPhasingDictionary(phased_df, tx_df)
  
  # Step 8: Create DoCo-level count
  cat("Step 8/8: Creating DoCo-level count...\n")
  doco_count <- blessy.createDoCoCount(phasing_dict, tx_count)
  
  # Step 9: Return the results
  cat("Pipeline completed. Returning results...\n")
  return(list(
    phasing_dict = phasing_dict,
    doco_count = doco_count
  ))
}
