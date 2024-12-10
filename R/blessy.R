#' blessy: A Wrapper Function to Run All Steps in the blessy's Pipeline
#'
#' This wrapper function performs the entire workflow of fetching domain and transcript tracks,
#' mapping domains to transcripts, deduplication, phasing, creating a phasing dictionary, 
#' and generating DoCo-level count data frames given a transcript count.
#'
#' @param genomeAssembly A string specifying the genome/assembly UCSC identifier (e.g., "hg38").
#' @param transcriptAnnotation A string or data frame specifying the transcript annotation track UCSC identifier (e.g., "wgEncodeGencodeBasicV44").
#'   If a data frame is provided, refer to `blessy.usingCustomAnnotation`.
#' @param domainAnnotation A string or data frame specifying the domain annotation track.
#'   If a data frame is provided, refer to `blessy.usingCustomAnnotation`.
#' @param transcriptCount A data frame containing transcript counts with:
#'   - \code{TranscriptID}: Transcript IDs (required as the first column).
#'   - Additional columns representing RNA-seq counts for biological samples (must be numeric).
#' @param coordinates Logical flag indicating whether to include genomic coordinates and strand information
#'   in the \code{DoCo} string during the phasing step. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#'   - \code{phasing_dict}: The phasing dictionary data frame containing information of Gene, DoCo, and Transcript.
#'   - \code{doco_count}: The DoCo-level count data frame created by aggregating transcripts of the same DoCo.
#'
#' @examples
#' # Example data for blessy
#' genomeAssembly <- "hg38"
#' 
#' transcriptAnnotation <- "wgEncodeGencodeBasicV44"
#' 
#' domainAnnotation <- "unipDomain"
#'
#' transcriptCount <- data.frame(
#'   TranscriptID = c("tx1", "tx2", "tx3"),
#'   Sample1 = c(10, 20, 30),
#'   Sample2 = c(15, 25, 35),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Run blessy with coordinates (default)
#' results_with_coords <- blessy(
#'   genomeAssembly, 
#'   transcriptAnnotation, 
#'   domainAnnotation, 
#'   transcriptCount
#' )
#'
#' # Run blessy without coordinates
#' results_no_coords <- blessy(
#'   genomeAssembly, 
#'   transcriptAnnotation, 
#'   domainAnnotation, 
#'   transcriptCount,
#'   coordinates = FALSE
#' )
#'
#' # Access the results
#' phasing_dict <- results_with_coords$phasing_dict
#' doco_count <- results_with_coords$doco_count
#'
#' @export
blessy <- function(genomeAssembly, transcriptAnnotation, domainAnnotation, transcriptCount, coordinates = TRUE) {
  # Step 1: Fetch transcript and domain annotation tracks
  cat("Step 1/9: Fetching transcript and domain annotation tracks...\n")
  tx_df <- blessy.getTranscriptTrack(genomeAssembly, transcriptAnnotation)
  domain_df <- blessy.getDomainTrack(genomeAssembly, domainAnnotation)
  
  # Step 2: Convert data frames to GRangesList objects
  cat("Step 2/9: Converting data frames to GRangesList objects...\n")
  tx_grangesList <- blessy.dfToGRangesList(tx_df)
  domain_grangesList <- blessy.dfToGRangesList(domain_df)
  
  # Step 3: Map domains to transcripts
  cat("Step 3/9: Mapping domains to transcripts...\n")
  mapped_df <- blessy.mapDomainToTranscript(tx_grangesList, domain_grangesList, tx_df, domain_df)
  
  # Step 4: Add exon and block starts/ends
  cat("Step 4/9: Adding exon and block starts/ends...\n")
  cat("Note: Patience is bitter, but its fruit is sweet. \n")
  starts_ends_df <- blessy.addStartsEnds(mapped_df)
  
  # Step 5: Deduplicate domain mappings
  cat("Step 5/9: Deduplicating domain mappings...\n")
  deduplicated_df <- blessy.domainDeduplication(starts_ends_df)
  
  # Step 6: Create phasing information with the 'coordinates' parameter
  cat(sprintf("Step 6/9: Creating phasing information (coordinates = %s)...\n", coordinates))
  phased_df <- blessy.domainPhasing(deduplicated_df, coordinates = coordinates)
  
  # Step 7: Create the phasing dictionary
  cat("Step 7/9: Creating the phasing dictionary...\n")
  phasing_dict <- blessy.createPhasingDictionary(phased_df, tx_df)
  
  # Step 8: Create DoCo-level count
  cat("Step 8/9: Creating DoCo-level count...\n")
  doco_count <- blessy.createDoCoCount(phasing_dict, transcriptCount)
  
  # Step 9: Return the results
  cat("Pipeline completed. Returning results...\n")
  return(list(
    phasing_dict = phasing_dict,
    doco_count = doco_count
  ))
}
