#' Convert a GTF Transcript annotation into a BED-like Data Frame 
#'
#' This function imports a GTF file containing transcript annotation using \pkg{rtracklayer} and re-formats
#' the data into a BED-like data frame, grouping by \code{transcript_id}.
#' Exons for each transcript are collapsed into blocks, producing a single row
#' per transcript. The resulting data frame is similar to a BED12 format.
#'
#' @param gtf_file A character string specifying the path to the GTF file.
#' @param feature_type A character string specifying which feature type to use
#'   when determining blocks (default "exon"). Typically, GTFs label exons with
#'   \code{"exon"}.
#' @param geneID A logical indicating whether the \code{geneName} column should
#'   take the value of \code{gene_id} (when \code{TRUE}) instead of \code{gene_name}.
#'   Default is \code{FALSE}, meaning \code{geneName} will be \code{gene_name} if
#'   available.
#'
#' @return A data frame (tibble) in BED-like format with the following columns:
#' \itemize{
#'   \item \strong{chrom} - Chromosome or scaffold name (from \code{seqnames}).
#'   \item \strong{chromStart} - Transcript start position (minimum of all exon starts).
#'   \item \strong{chromEnd} - Transcript end position (maximum of all exon ends).
#'   \item \strong{name} - The \code{transcript_id}.
#'   \item \strong{score} - A numeric score (set to 0 by default).
#'   \item \strong{strand} - The strand (\code{+} or \code{-}).
#'   \item \strong{thickStart} - Currently set to the same as \code{chromStart}.
#'   \item \strong{thickEnd} - Currently set to the same as \code{chromEnd}.
#'   \item \strong{itemRgb} - Color value, default \code{"0,0,0"}.
#'   \item \strong{blockCount} - Number of exons (blocks).
#'   \item \strong{blockSizes} - A comma-separated string of exon lengths.
#'   \item \strong{blockStarts} - A comma-separated string of exon start offsets
#'         relative to \code{chromStart}.
#'   \item \strong{geneName} - Either \code{gene_name} (default) or \code{gene_id}
#'         (if \code{geneID} = \code{TRUE}).
#' }
#'
#' @details
#' This is useful when you want a compact representation of transcripts similar
#' to BED12 format. Each transcript occupies a single row, and its individual
#' exons are shown as "blocks".
#'
#' @import rtracklayer
#' @importFrom dplyr filter group_by arrange summarize first mutate select
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert a GTF file to a BED-like data frame:
#' bed_df <- blessy.GTFtoDataFrame("path/to/gencode.v38.annotation.gtf", geneID = FALSE)
#' head(bed_df)
#' }
blessy.GTFtoBEDdf <- function(gtf_file, feature_type = "exon", geneID = FALSE) {
  
  # 1) Load the GTF via rtracklayer
  gtf_gr <- rtracklayer::import(gtf_file)
  gtf_df <- as.data.frame(gtf_gr)
  
  # 2) Optionally filter to only the desired feature type (usually "exon")
  library(dplyr)
  gtf_df <- gtf_df %>%
    filter(type == feature_type)
  
  # 3) Group by transcript_id, then collect exons into BED-like blocks
  bed_df <- gtf_df %>%
    # Keep only rows with a valid transcript_id
    filter(!is.na(transcript_id)) %>%
    group_by(seqnames, transcript_id) %>%
    arrange(start, .by_group = TRUE) %>%
    summarize(
      chrom       = first(seqnames),
      chromStart  = min(start),
      chromEnd    = max(end),
      name        = first(transcript_id),
      score       = 0,  # No score info in typical GTF; set to 0
      strand      = first(strand),
      thickStart  = min(start),
      thickEnd    = max(end),
      itemRgb     = "0,0,0",
      blockCount  = n(),  # number of exons for this transcript
      blockSizes  = paste(end - start, collapse = ","),
      blockStarts = paste(start - min(start), collapse = ","),
      geneName    = if (geneID) {
        # If geneID = TRUE, pull from gene_id
        if ("gene_id" %in% names(.)) first(gene_id) else NA_character_
      } else {
        # Otherwise, pull from gene_name if present, else NA
        if ("gene_name" %in% names(.)) first(gene_name) else NA_character_
      },
      .groups     = "drop"
    )
  
  # 4) Finalize column order to mimic a BED-like structure
  bed_df <- bed_df %>%
    select(
      chrom, chromStart, chromEnd, name, score, strand,
      thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts, geneName
    )
  
  # 5) Return the result
  return(bed_df)
}
