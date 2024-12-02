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
#'   - \code{chromStarts}: Comma-separated relative block start positions.
#'   - \code{blockSizes}: Comma-separated block sizes.
#'
#' @return The input data frame with four new columns appended:
#'   - \code{exonStarts}: Computed exon start positions (comma-separated).
#'   - \code{exonEnds}: Computed exon end positions (comma-separated).
#'   - \code{blockStarts}: Computed block start positions (comma-separated).
#'   - \code{blockEnds}: Computed block end positions (comma-separated).
#'
#' @importFrom dplyr select mutate group_by summarise left_join %>%
#' @importFrom tidyr unnest_longer
#'
#' @examples
#' # Example data frame
#' df <- data.frame(
#'   txStart = c(1000, 2000),
#'   exonRelativeStarts = c("0,100,200", "0,50,100"),
#'   exonSizes = c("100,50,25", "80,40,30"),
#'   chromStart = c(1000, 2000),
#'   chromStarts = c("0,150,300", "0,100,200"),
#'   blockSizes = c("100,50,25", "80,40,30")
#' )
#'
#' # Apply the function to compute start and end positions
#' df <- blessy.addStartsEnds(df)
#'
#' @export
blessy.addStartsEnds <- function(df) {
  # Ensure required columns
  required_columns <- c("txStart", "exonRelativeStarts", "exonSizes", "chromStart", "chromStarts", "blockSizes")
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing:", paste(missing_columns, collapse = ", ")))
  }
  
  # Add a unique identifier for each row
  df$row_id <- seq_len(nrow(df))
  
  # Process exonStarts and exonEnds
  exon_df <- df %>%
    select(row_id, txStart, exonRelativeStarts, exonSizes) %>%
    mutate(
      exonRelativeStarts = strsplit(as.character(exonRelativeStarts), ","),
      exonSizes = strsplit(as.character(exonSizes), ",")
    ) %>%
    unnest_longer(c(exonRelativeStarts, exonSizes)) %>%
    mutate(
      exonStarts = as.numeric(exonRelativeStarts) + txStart,
      exonEnds = exonStarts + as.numeric(exonSizes)
    ) %>%
    group_by(row_id) %>%
    summarise(
      exonStarts = paste(exonStarts, collapse = ","),
      exonEnds = paste(exonEnds, collapse = ","),
      .groups = "drop"
    )
  
  # Process blockStarts and blockEnds
  block_df <- df %>%
    select(row_id, chromStart, chromStarts, blockSizes) %>%
    mutate(
      chromStarts = strsplit(as.character(chromStarts), ","),
      blockSizes = strsplit(as.character(blockSizes), ",")
    ) %>%
    unnest_longer(c(chromStarts, blockSizes)) %>%
    mutate(
      blockStarts = as.numeric(chromStarts) + chromStart,
      blockEnds = blockStarts + as.numeric(blockSizes)
    ) %>%
    group_by(row_id) %>%
    summarise(
      blockStarts = paste(blockStarts, collapse = ","),
      blockEnds = paste(blockEnds, collapse = ","),
      .groups = "drop"
    )
  
  # Merge the computed columns back into the original data frame
  df <- df %>%
    left_join(exon_df, by = "row_id") %>%
    left_join(block_df, by = "row_id") %>%
    select(-row_id)  # Remove the temporary row identifier
  
  return(df)
}
