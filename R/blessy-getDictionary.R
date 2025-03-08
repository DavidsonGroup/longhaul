#' blessy.getDictionary: A Wrapper Function to Run All Steps in the blessy's Pipeline related to dictionary creation
#'
#' This wrapper function performs the workflow of fetching domain and transcript tracks,
#' mapping domains to transcripts, deduplication, phasing, creating a phasing dictionary.
#'
#' @param genomeAssembly A string specifying the genome/assembly UCSC identifier (e.g., "hg38").
#' @param transcriptAnnotation A string or data frame specifying the transcript annotation track UCSC identifier (e.g., "wgEncodeGencodeBasicV44").
#'   If a data frame is provided, refer to `blessy.usingCustomAnnotation`.
#' @param domainAnnotation A string or data frame specifying the domain annotation track.
#'   If a data frame is provided, refer to `blessy.usingCustomAnnotation`.
#' @param unique_domain Logical flag indicating whether to deduplicate domains by keeping only unique domains per transcript. 
#'   Defaults to \code{FALSE}.
#' @param coordinates Logical flag indicating whether to include genomic coordinates and strand information
#'   in the \code{DoCo} string during the phasing step. Defaults to \code{TRUE}.
#'
#' @return A data frame
#'   - \code{phasing_dict}: The phasing dictionary data frame containing information of Gene, DoCo, and Transcript.
#'
#' @examples
#' # Example data for blessy
#' genomeAssembly <- "hg38"
#' 
#' transcriptAnnotation <- "wgEncodeGencodeBasicV44"
#' 
#' domainAnnotation <- "unipDomain"
#'
#'
#' # Run blessy with coordinates (default)
#' dict_with_coords <- blessy.getDictionary(
#'   genomeAssembly, 
#'   transcriptAnnotation, 
#'   domainAnnotation, 
#' )
#'
#' # Run blessy without coordinates
#' dict_no_coords <- blessy.getDictionary(
#'   genomeAssembly, 
#'   transcriptAnnotation, 
#'   domainAnnotation, 
#'   coordinates = FALSE
#' )
#'
#'
#' @export
blessy.getDictionary <- function(genomeAssembly, transcriptAnnotation, domainAnnotation, unique_domain = FALSE, coordinates = TRUE) {
  # Step 1: Fetch transcript and domain annotation tracks
  cat("Step 1/6: Fetching transcript and domain annotation tracks...\n")
  tx_df <- blessy.getTranscriptTrack(genomeAssembly, transcriptAnnotation)
  domain_df <- blessy.getDomainTrack(genomeAssembly, domainAnnotation)
  
  # Step 2: Map domains to transcripts
  cat("Step 2/6: Matching domains to transcripts...\n")
  mapped_df <- blessy.mapDomainToTranscript(tx_df, domain_df)
  
  # Step 3: Add exon and block starts/ends
  cat("Step 3/6: Adding exon and block starts/ends...\n")
  cat("Note: Patience is bitter, but its fruit is sweet. \n")
  starts_ends_df <- blessy.addStartsEnds(mapped_df)
  
  # Step 4: Deduplicate domain mappings
  cat("Step 4/6: Deduplicating domain mappings...\n")
  deduplicated_df <- blessy.domainDeduplication(starts_ends_df, unique_domain = unique_domain)
  
  # Step 5: Create phasing information with the 'coordinates' parameter
  cat("Step 5/6: Creating phasing information...\n")
  phased_df <- blessy.domainPhasing(deduplicated_df, coordinates = coordinates)
  
  # Step 6: Create the phasing dictionary
  cat("Step 6/6: Creating the phasing dictionary...\n")
  phasing_dict <- blessy.createPhasingDictionary(phased_df, tx_df)
  
  # Step 7: Return the results
  cat("Dictionary creation pipeline completed. Returning results...\n")
  return(phasing_dict)
 
}
