#' Phase Domains for Transcript-Domain Mapping Data Frame
#'
#' This function processes a transcript-domain data frame, phasing domains along transcripts
#' by sorting based on strand and location to create a concatenated domain-combination (`DoCo`) string.
#'
#' @param tx_domain_df A data frame containing transcript-domain data with the following required columns:
#'   - \code{Transcript}: Transcript identifier.
#'   - \code{Domain}: Domain identifier (may include NA values).
#'   - \code{chrom}: Chromosome identifier for the domain.
#'   - \code{chromStart}, \code{chromEnd}: Start and end positions for the domain.
#'   - \code{strand}: Strand information for the domain.
#'   - \code{Gene}: Gene identifier.
#'
#' @return The input data frame with an additional column, \code{DoCo}, which contains the domain-combination
#'  string for each transcript.
#'
#' @importFrom dplyr mutate group_by arrange ungroup first
#' @importFrom stats na.omit
#'
#' @examples
#' # Example input data frame
#' mapping_df <- data.frame(
#'   Transcript = c("tx1", "tx1", "tx2"),
#'   Domain = c("domainA", "domainB", NA),
#'   chrom2 = c("chr1", "chr1", "chr2"),
#'   chromStart = c(1000, 2000, 3000),
#'   chromEnd = c(1100, 2100, 3100),
#'   strand2 = c("+", "-", "+"),
#'   strand = c("+", "+", "-"),
#'   Gene = c("geneA", "geneA", "geneB")
#' )
#'
#' # Apply domain phasing
#' phased_df <- blessy.domainPhasing(mapping_df)
#'
#' @export
blessy.domainPhasing <- function(tx_domain_df) {
  # Ensure required columns exist
  required_columns <- c("Transcript", "Domain", "chrom", "chromStart", "chromEnd", "strand", "Gene")
  missing_columns <- setdiff(required_columns, colnames(tx_domain_df))
  if (length(missing_columns) > 0) {
    stop("The dataframe must contain the following columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # Ensure chromStart and chromEnd are numeric for proper sorting
  tx_domain_df <- tx_domain_df %>%
    mutate(
      chromStart = as.numeric(chromStart),
      chromEnd = as.numeric(chromEnd)
    )
  
  # Group and process data
  tx_domain_df <- tx_domain_df %>%
    group_by(Transcript) %>%
    arrange(
      if (first(strand) == "+") chromStart else -chromStart,
      .by_group = TRUE
    ) %>%
    mutate(
      DoCo = if (all(is.na(Domain))) {
        paste0(";;; ", first(Gene))
      } else {
        # Create domain info strings
        domain_info <- ifelse(
          is.na(Domain),
          NA_character_,
          paste0(
            Domain, "::", chrom, ":", chromStart, "-", chromEnd, "(", strand, ")"
          )
        )
        # Concatenate domain info strings
        domain_string <- paste(na.omit(domain_info), collapse = ",")
        paste0(domain_string, ";;; ", first(Gene))
      }
    ) %>%
    ungroup()
  
  return(tx_domain_df)
}
