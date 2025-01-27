#' Create DoCo-Level Count 
#'
#' This function aggregates RNA-seq counts from a count data frame at the DoCo level. 
#' Transcripts are grouped based on their DoCo assignment from a dictionary data frame (dict),
#' and counts are summed across biological samples. Unmatched transcripts are grouped into the DoCo
#' class ";;;". 
#'
#' @param dict A data frame containing the Transcript and DoCo columns:
#'   - \code{Transcript}: Transcript identifier.
#'   - \code{DoCo}: Domain-combination string (grouping information).
#' @param count_df A data frame containing RNA-seq counts with:
#'   - \code{rownames}: Transcript IDs.
#'   - Columns representing RNA-seq counts for biological samples (must be numeric).
#'
#' @return A data frame aggregated at the DoCo level with the following columns:
#'   - \code{DoCo}: Domain-combination string.
#'   - RNA-seq counts for biological samples.
#'
#' @importFrom dplyr mutate
#' @importFrom stats aggregate
#'
#' @examples
#' # Example dictionary data frame
#' dict <- data.frame(
#'   Transcript = c("Tx1", "Tx2", "Tx3"),
#'   DoCo = c("D1,D2;;; GeneA", "D3;;; GeneB", ";;; GeneC"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Example count data frame
#' count_df <- data.frame(
#'   Sample1 = c(10, 20, 5, 7),
#'   Sample2 = c(15, 25, 8, 10),
#'   row.names = c("Tx1", "Tx2", "Tx4", "Tx5")
#' )
#'
#' # Create DoCo-level counts
#' doco_count <- blessy.createDoCoCount(dict, count_df)
#'
#' @export
blessy.createDoCoCount <- function(dict, count_df) {
  # Check for required columns in dict
  required_dict_columns <- c("Transcript", "DoCo")
  if (!all(required_dict_columns %in% colnames(dict))) {
    stop("The dictionary dataframe must contain the following columns: ", paste(required_dict_columns, collapse = ", "))
  }
  
  # Ensure rownames are present in count_df
  if (is.null(rownames(count_df))) {
    stop("The count_df must have rownames representing Transcript IDs.")
  }
  
  # Check that all columns in count_df are numeric
  if (!all(sapply(count_df, is.numeric))) {
    stop("Non-numeric values detected in count columns.")
  }
  
  # Remove version numbers from Transcript and rownames
  if (any(grepl("\\.[0-9]+$", dict$Transcript))) {
    warning("Version numbers found in dict$Transcript. Removing them...")
    dict$Transcript <- sub("\\.[0-9]+$", "", dict$Transcript)
  }
  
  if (any(grepl("\\.[0-9]+$", rownames(count_df)))) {
    warning("Version numbers found in rownames of count_df. Removing them...")
    rownames(count_df) <- sub("\\.[0-9]+$", "", rownames(count_df))
  }
  
  # Convert rownames of count_df into a column for merging
  count_df$TranscriptID <- rownames(count_df)
  rownames(count_df) <- NULL
  
  # Merge dict with count_df to map transcripts to DoCo
  merged_df <- merge(count_df, dict, by.x = "TranscriptID", by.y = "Transcript", all.x = TRUE)
  
  # Handle unmatched transcripts
  unmatched <- is.na(merged_df$DoCo)
  num_unmatched <- sum(unmatched)
  
  if (num_unmatched > 0) {
    # Assign unmatched transcripts to new DoCo ";;;"
    merged_df$DoCo[unmatched] <- ";;;"
    warning("Uncommon transcripts found, grouping them into the empty ';;;' DoCo group.")
    warning(paste("Number of uncommon transcripts grouped:", num_unmatched))
  }
  
  # Aggregate counts at the DoCo level
  sample_columns <- colnames(count_df)[-ncol(count_df)]
  aggregated_df <- aggregate(
    merged_df[, sample_columns],
    by = list(DoCo = merged_df$DoCo),
    FUN = sum
  )
  
  # Ensure the first column (DoCo) is set as row names
  rownames(aggregated_df) <- aggregated_df$DoCo
  
  # Drop the 'DoCo' column 
  aggregated_df <- aggregated_df[, -1]
  
  # Return the modified data frame
  return(aggregated_df)
}