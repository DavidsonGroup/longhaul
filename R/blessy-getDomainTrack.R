#' Fetch and Convert UCSC Domain-type Track to Data Frame
#'
#' This function fetches a domain annotation track from the UCSC Genome Browser using UCSC.utils
#' and converts it into a BED-formatted data frame in R, notably the UniProt and Pfam tracks.
#'
#' @param genome A string specifying the genome/assembly identifier (e.g., "hg38").
#' @param track A string specifying the table identifier of the domain track (e.g., "unipDomain").
#'
#' @return A data frame in BED format containing the domain annotation track.
#'
#' @importFrom UCSC.utils fetch_UCSC_track_data
#' @importFrom dplyr %>%
#'
#' @examples
#' # Example: Fetch the "unipDomain" annotation track
#' bed_df <- blessy.getDomainTrack("hg38", "unipDomain")
#'
#' @export
blessy.getDomainTrack <- function(genome, track) {
  # Step 1: Fetch the UCSC track data using UCSC.utils
  track_data <- fetch_UCSC_track_data(genome, track)
  
  # Step 2: Rename necessary columns to match the BED format
  colnames(track_data) <- colnames(track_data) %>%
    gsub("reserved", "itemRgb", .) %>%      # Rename reserved to itemRgb
    gsub("chromStarts", "blockStarts", .)  # Rename chromStarts to blockStarts
  
  # Step 3: Define the required 12 columns for a valid BED file
  required_columns <- c("chrom", "chromStart", "chromEnd", "name", "score",
                        "strand", "thickStart", "thickEnd", "itemRgb",
                        "blockCount", "blockSizes", "blockStarts")
  
  # Step 4: Check if any required columns are missing, and fill with defaults if necessary
  for (col in required_columns) {
    if (!col %in% colnames(track_data)) {
      track_data[[col]] <- if (col == "itemRgb") "0,0,0" else NA  # Default for itemRgb, NA otherwise
    }
  }
  
  # Step 5: Set thickStart and thickEnd to chromStart and chromEnd respectively
  track_data$thickStart <- track_data$chromStart
  track_data$thickEnd <- track_data$chromEnd
  
  # Step 6: Select and reorder only the required columns
  bed_df <- track_data[, required_columns, drop = FALSE]
  
  # Step 7: Return the BED-formatted data frame
  return(bed_df)
}
