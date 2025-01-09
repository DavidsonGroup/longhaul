#' Convert a BED-like Data Frame to a CompressedGRangesList
#'
#' This function converts a BED-like data frame into a `CompressedGRangesList` object, where each row of the initial 
#' data frame is represented as a `GRanges` object in the returned object.
#'
#' @param df A data frame containing BED-like annotations. Required columns include:
#'   - \code{chrom}: Chromosome names.
#'   - \code{thickStart}: Start positions of genomic ranges (numeric).
#'   - \code{thickEnd}: End positions of genomic ranges (numeric).
#'   - \code{strand}: Strand information (+, -, or *).
#'   Optional additional metadata columns can be included.
#'
#' @return A `CompressedGRangesList` object where each element corresponds to a row in the input data frame.
#'
#' @import GenomicRanges
#'
#' @examples
#' # Example BED-like data frame
#' bed_df <- data.frame(
#'   chrom = c("chr1", "chr1", "chr2"),
#'   thickStart = c(1000, 2000, 3000),
#'   thickEnd = c(1500, 2500, 3500),
#'   strand = c("+", "-", "+"),
#'   name = c("feature1", "feature2", "feature3")
#' )
#'
#' # Convert the data frame to a GRangesList
#' grl <- blessy.dfToGRanges(bed_df)
#' 
#' @export
blessy.dfToGRangesList <- function(df) {
  # Ensure required columns exist
  required_columns <- c("chrom", "thickStart", "thickEnd", "strand")
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(
      paste(
        "The following required columns are missing:",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  # Convert required columns to appropriate types
  df$thickStart <- as.numeric(df$thickStart)
  df$thickEnd <- as.numeric(df$thickEnd)
  df$strand <- as.character(df$strand)
  
  # Create a GRanges object directly from the data frame
  gr <- makeGRangesFromDataFrame(
    df,
    seqnames.field = "chrom",
    start.field = "thickStart",
    end.field = "thickEnd",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )
  
  # Split the GRanges object into a GRangesList where each element is a single GRanges
  grl <- split(gr, seq_along(gr))
  
  # Optionally set names for the GRangesList elements
  if (!is.null(rownames(df))) {
    names(grl) <- rownames(df)
  } else if ("id" %in% colnames(df)) {
    names(grl) <- df$id
  } else {
    names(grl) <- seq_along(gr)
  }
  
  return(grl)
}
