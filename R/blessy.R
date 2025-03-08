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
#' @param transcriptCount A data frame containing RNA-seq counts with:
#'   - \code{rownames}: Transcript IDs.
#'   - Columns representing RNA-seq counts for biological samples (must be numeric).
#' @param unique_domain Logical flag indicating whether to deduplicate domains by keeping only unique domains per transcript. 
#'   Defaults to \code{FALSE}.
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
#'   Sample1 = c(10, 20, 30),
#'   Sample2 = c(15, 25, 35),
#'   row.names = c("ENST00000000412", "ENST00000000442", "ENST00000001008")
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

blessy <- function(genomeAssembly, transcriptAnnotation,
       domainAnnotation, transcriptCount,
       unique_domain = FALSE, coordinates = TRUE) {

  
  cat("Step 1-6/7: Creating DoCo-dictionary...\n")
  phasing_dict<-blessy.getDictionary(genomeAssembly,transcriptAnnotation,
  			             domainAnnotation, unique_domain, coordinates)

  cat("Step 7/7: Creating DoCo-level count...\n")
  doco_count <- blessy.createDoCoCount(phasing_dict, transcriptCount)
  
  cat("Pipeline completed. Returning results...\n")
  return(list(
    phasing_dict = phasing_dict,
    doco_count = doco_count
  ))
}
