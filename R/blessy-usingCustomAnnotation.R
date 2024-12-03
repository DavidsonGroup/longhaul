#' blessy's Wrapper Function using Custom Annotation
#'
#' This wrapper function performs the workflow starting from mapping domains to transcripts
#' using user-provided transcript and domain annotation data frames of BED-like format. It creates a phasing dictionary
#' and generates DoCo-level count data frames.
#'
#' @param tx_df A BED-like data frame representing transcript annotation. The annotation must include information on exon locations in `blockCount`, `blockSizes`, `blockStarts` and gene name in `geneName`.
#' @param domain_df A BED-like data frame representing domain annotation. The annotation must include information on domain block locations in in `blockCount`, `blockSizes`, `blockStarts`.
#' @param tx_count A data frame containing transcript counts with:
#'   - \code{TranscriptID}: Transcript IDs (required as the first column).
#'   - Additional columns representing RNA-seq counts for biological samples (must be numeric).
#'
#' @return A list containing:
#'   - \code{phasing_dict}: The phasing dictionary data frame
#'   - \code{doco_count}: The DoCo-level count data frame
#'
#'
#' @export
blessy.usingCustomAnnotation <- function(tx_df, domain_df, tx_count) {
  #Step 2: Convert dataframes to GRanges object
  tx_grangesList <- blessy.dfToGRangesList(tx_df)
  domain_grangesList <- blessy.dfToGRangesList(domain_df)
  
  # Step 3: Map domains to transcripts
  mapped_df <- blessy.mapDomainToTranscript(tx_grangesList, domain_grangesList, tx_df, domain_df)
  
  # Step 4: Add exon and block starts/ends
  starts_ends_df <- blessy.addStartsEnds(mapped_df)
  
  # Step 5: Deduplicate domain mappings
  deduplicated_df <- blessy.domainDeduplication(starts_ends_df)
  
  # Step 6: Create phasing information
  phased_df <- blessy.domainPhasing(deduplicated_df)
  
  # Step 7: Create the phasing dictionary
  phasing_dict <- blessy.createPhasingDictionary(phased_df, tx_df)
  
  # Step 8: Create DoCo-level counts
  doco_count <- blessy.createDoCoCount(phasing_dict, tx_count)
  
  # Return the phasing dictionary and DoCo counts as a list
  return(list(
    phasing_dict = phasing_dict,
    doco_count = doco_count
  ))
}
