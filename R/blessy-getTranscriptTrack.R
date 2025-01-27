#' Fetch and Convert UCSC Transcript-type Track to Data Frame
#'
#' This function fetches a transcript annotation track from the UCSC Genome Browser using UCSC.utils and
#' converts it into a BED-formatted data frame in R, notably the NCBI RefSeq or GENCODE tracks.
#'
#' @param genome A string specifying the genome/assembly identifier (e.g., "hg38").
#' @param track A string specifying the table identifier of the annotation track (e.g., "wgEncodeGencodeBasicV44").
#'
#' @return A data frame in BED format containing the transcript annotation track.
#'
#' @import UCSC.utils
#'
#' @examples
#' # Example: Fetch the "wgEncodeGencodeBasicV44" annotation transcript annotation track
#' bed_df <- blessy.getTranscriptTrack("hg38", "wgEncodeGencodeBasicV44")
#'
#' @export
blessy.getTranscriptTrack <- function(genome, track) {
  # Step 1: Fetch the UCSC track data using UCSC.utils
  track_data <- fetch_UCSC_track_data(genome, track)
  
  # Step 2: Ensure necessary columns are present in the fetched data
  required_columns <- c("chrom", "txStart", "txEnd", "name", "score", "strand",
                        "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "name2")
  
  if (!all(required_columns %in% names(track_data))) {
    stop("The fetched data is missing one or more required columns.")
  }
  
  # Step 3: Convert the track data to BED format
  bed_df <- data.frame(
    chrom = track_data$chrom,
    chromStart = track_data$txStart,
    chromEnd = track_data$txEnd,
    name = track_data$name,
    score = ifelse(is.na(track_data$score), 0, track_data$score), # Fill in missing scores with 0
    strand = track_data$strand,
    thickStart = track_data$cdsStart,
    thickEnd = track_data$cdsEnd,
    itemRgb = "0,0,0", # Use default RGB value
    blockCount = track_data$exonCount,
    blockSizes = apply(track_data, 1, function(row) {
      # Calculate block sizes by subtracting exonStarts from exonEnds
      starts <- as.numeric(unlist(strsplit(row["exonStarts"], ",")))
      ends <- as.numeric(unlist(strsplit(row["exonEnds"], ",")))
      paste(ends - starts, collapse = ",")
    }),
    blockStarts = apply(track_data, 1, function(row) {
      # Calculate block starts relative to txStart
      starts <- as.numeric(unlist(strsplit(row["exonStarts"], ",")))
      paste(starts - as.numeric(row["txStart"]), collapse = ",")
    }),
    geneName = track_data$name2 # Extra column
  )
  
  # Step 4: Return the BED-formatted data frame
  return(bed_df)
}