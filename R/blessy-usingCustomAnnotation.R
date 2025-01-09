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
#' @param coordinates Logical flag indicating whether to include genomic coordinates and strand information
#'   in the \code{DoCo} string during the phasing step. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#'   - \code{phasing_dict}: The phasing dictionary data frame.
#'   - \code{doco_count}: The DoCo-level count data frame.
#'
#' @export
blessy.usingCustomAnnotation <- function(tx_df, domain_df, tx_count, unique_domain = FALSE, coordinates = TRUE) {
  # Step 2: Check input types and convert to GRangesList if necessary
  if (is(tx_df, "GRangesList") && is(domain_df, "GRangesList")) {
    cat("Step 2/9: Input parameters are already GRangesList objects. Skipping conversion...\n")
    tx_grangesList <- tx_df
    domain_grangesList <- domain_df
  } else {
    cat("Step 2/9: Converting data frames to GRangesList objects...\n")
    tx_grangesList <- blessy.dfToGRangesList(tx_df)
    domain_grangesList <- blessy.dfToGRangesList(domain_df)
  }
  
  # Step 3: Map domains to transcripts
  cat("Step 3/9: Mapping domains to transcripts...\n")
  mapped_df <- blessy.mapDomainToTranscript(tx_grangesList, domain_grangesList, tx_df, domain_df)
  
  # Step 4: Add exon and block starts/ends
  cat("Step 4/9: Adding exon and block starts/ends...\n")
  cat("Note: Patience is bitter, but its fruit is sweet. \n")
  starts_ends_df <- blessy.addStartsEnds(mapped_df)
  
  # Step 5: Deduplicate domain mappings
  cat("Step 5/9: Deduplicating domain mappings...\n")
  deduplicated_df <- blessy.domainDeduplication(starts_ends_df, unique_domain = unique_domain)
  
  # Step 6: Create phasing information with the 'coordinates' parameter
  cat(sprintf("Step 6/9: Creating phasing information (coordinates = %s)...\n", coordinates))
  phased_df <- blessy.domainPhasing(deduplicated_df, coordinates = coordinates)
  
  # Step 7: Create the phasing dictionary
  cat("Step 7/9: Creating the phasing dictionary...\n")
  phasing_dict <- blessy.createPhasingDictionary(phased_df, tx_df)
  
  # Step 8: Create DoCo-level count
  cat("Step 8/9: Creating DoCo-level count...\n")
  doco_count <- blessy.createDoCoCount(phasing_dict, tx_count)
  
  # Return the results
  cat("Pipeline completed. Returning results...\n")
  return(list(
    phasing_dict = phasing_dict,
    doco_count = doco_count
  ))
}
