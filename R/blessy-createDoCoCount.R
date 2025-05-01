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
#' @importFrom data.table as.data.table .SD
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
  
  # Check numeric content
  if (is.data.frame(count_df)) {
     if (!all(sapply(count_df, is.numeric))) {
     stop("Non-numeric values detected in count_df.")
  }
  } else if (inherits(count_df, "Matrix")) {
    # OK: sparse matrix from Matrix package
    } else if (is.matrix(count_df)) {
    if (!is.numeric(count_df)) {
       stop("Base matrix is not numeric.")
  }
  } else {
    stop("count_df must be a data.frame, matrix, or sparse Matrix.")
    }

  # Remove transcript version numbers
  dict$Transcript <- sub("\\.[0-9]+$", "", dict$Transcript)
  rownames(count_df) <- sub("\\.[0-9]+$", "", rownames(count_df))
  if (any(grepl("\\.[0-9]+$", rownames(count_df)))) {
    warning("Version numbers found in rownames of count_df. Removing them...")
    rownames(count_df) <- sub("\\.[0-9]+$", "", rownames(count_df))
  }

  # Convert count_df to data.table with transcript IDs
  transcript_ids <- rownames(count_df)

  # Merge DoCo annotations
  matched <- match(transcript_ids, dict$Transcript)
  doco_vec <- dict$DoCo[matched]
  
  # Handle unmatched transcripts
  num_unmatched <- sum(is.na(doco_vec))
  if (num_unmatched > 0) {
    warning("Uncommon transcripts found, grouping them into the empty ';;;' DoCo group.")
    warning(paste("Number of uncommon transcripts grouped:", num_unmatched))
    doco_vec[is.na(doco_vec)] <- ";;;"
    }

    # Convert doco_vec to factor for grouping
    doco_factor <- factor(doco_vec)


    if (inherits(count_df, "Matrix")) {
       group_mat <- Matrix::sparse.model.matrix(~ doco + 0, data = data.frame(doco = doco_factor))
       aggregated_mat <- t(group_mat) %*% count_df
       rownames(aggregated_mat) <- levels(doco_factor)
       return(aggregated_mat)
    } #otherwise it's a regular data.frame      
    aggregated_df <- rowsum(count_df, group = doco_factor)
    return(aggregated_df)  # return as data.frame
}
