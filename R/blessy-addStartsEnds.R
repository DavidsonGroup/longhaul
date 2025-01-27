#' Add Exon and Block Start/End Positions to a Data Frame
#'
#' This function computes exon and block's start and end positions based on
#' BED-like annotation columns in the input data frame.
#' It returns the data frame with the computed columns appended.
#' 
#'
#' @param df A BED-like data frame containing required columns:
#'   - \code{txStart}: Transcript start position.
#'   - \code{exonRelativeStarts}: Comma-separated relative exon start positions.
#'   - \code{exonSizes}: Comma-separated exon sizes.
#'   - \code{chromStart}: Chromosomal start position.
#'   - \code{domainBlockStarts}: Comma-separated relative block start positions.
#'   - \code{blockSizes}: Comma-separated block sizes.
#'
#' @return The input data frame with four new columns appended:
#'   - \code{exonStarts}: Computed exon start positions (comma-separated).
#'   - \code{exonEnds}: Computed exon end positions (comma-separated).
#'   - \code{blockStarts}: Computed block start positions (comma-separated).
#'   - \code{blockEnds}: Computed block end positions (comma-separated).
#'
#'
#' @examples
#' # Example data frame
#' df <- data.frame(
#'   txStart = c(1000, 2000),
#'   exonRelativeStarts = c("0,100,200", "0,50,100"),
#'   exonSizes = c("100,50,25", "80,40,30"),
#'   chromStart = c(1000, 2000),
#'   domainBlockStarts = c("0,150,300", "0,100,200"),
#'   blockSizes = c("100,50,25", "80,40,30")
#' )
#'
#' # Apply the function to compute start and end positions
#' df <- blessy.addStartsEnds(df)
#'
#' @export
blessy.addStartsEnds <- function(df) {
  # Ensure required columns
  required_columns <- c("txStart", "exonRelativeStarts", "exonSizes", 
                        "chromStart", "domainBlockStarts", "blockSizes")
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing:", paste(missing_columns, collapse = ", ")))
  }
  
  # Split all exon and block data up front
  exon_rel_starts_list <- strsplit(as.character(df$exonRelativeStarts), ",")
  exon_sizes_list <- strsplit(as.character(df$exonSizes), ",")
  
  block_starts_list <- strsplit(as.character(df$domainBlockStarts), ",")
  block_sizes_list <- strsplit(as.character(df$blockSizes), ",")
  
  # Compute exon starts/ends
  # mapply will iterate over each row, allowing vectorized arithmetic inside the function:
  exon_results <- mapply(function(tx_start, rel_starts, sizes) {
    rel_starts_num <- as.numeric(rel_starts)
    sizes_num <- as.numeric(sizes)
    exon_starts <- rel_starts_num + tx_start
    exon_ends <- exon_starts + sizes_num
    list(
      exonStarts = paste(exon_starts, collapse = ","),
      exonEnds = paste(exon_ends, collapse = ",")
    )
  }, 
  tx_start = df$txStart, 
  rel_starts = exon_rel_starts_list, 
  sizes = exon_sizes_list, 
  SIMPLIFY = FALSE)
  
  # Compute block starts/ends
  block_results <- mapply(function(chrom_start, rel_starts, sizes) {
    rel_starts_num <- as.numeric(rel_starts)
    sizes_num <- as.numeric(sizes)
    block_starts <- rel_starts_num + chrom_start
    block_ends <- block_starts + sizes_num
    list(
      blockStarts = paste(block_starts, collapse = ","),
      blockEnds = paste(block_ends, collapse = ",")
    )
  },
  chrom_start = df$chromStart,
  rel_starts = block_starts_list,
  sizes = block_sizes_list,
  SIMPLIFY = FALSE)
  
  # Extract computed values into new columns
  df$exonStarts <- sapply(exon_results, `[[`, "exonStarts")
  df$exonEnds <- sapply(exon_results, `[[`, "exonEnds")
  
  df$blockStarts <- sapply(block_results, `[[`, "blockStarts")
  df$blockEnds <- sapply(block_results, `[[`, "blockEnds")
  
  return(df)
}
