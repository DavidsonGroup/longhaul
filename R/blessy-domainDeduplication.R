#' Domain Deduplication for Transcript-Domain Mapping
#'
#' This function deduplicates a transcript-domain mapping data frame by filtering out 
#' domains with non-exact matching. It ensures that only matches with relevant
#' exon-block relationships are retained.
#'
#' @param tx_domain_df A data frame containing transcript-domain mapping data, including columns:
#'   - \code{exonStarts}, \code{exonEnds}: Comma-separated exon start and end positions.
#'   - \code{blockStarts}, \code{blockEnds}: Comma-separated block start and end positions.
#'
#' @return A deduplicated data frame of transcript-domain mapping.
#'
#' @import dplyr 
#' @import tidyr
#'
#' @examples
#' # Example data frame
#' intersect_df <- data.frame(
#'   exonStarts = c("1000,2000,3000", "1500,2500,3500"),
#'   exonEnds = c("1100,2100,3100", "1600,2600,3600"),
#'   blockStarts = c("1000,2000,3000", "1500,2500,3500"),
#'   blockEnds = c("1100,2100,3100", "1600,2600,3600")
#' )
#'
#' # Apply the function to deduplicate domains
#' refined_df <- blessy.domainDeduplication(intersect_df)
#'
#' @export
blessy.domainDeduplication <- function(tx_domain_df) {
  # Add a unique identifier for each row/domain
  tx_domain_df <- tx_domain_df %>% mutate(row_id = row_number())
  
  # Split and unnest exonStarts and exonEnds into long format
  exon_df <- tx_domain_df %>%
    select(row_id, exonStarts, exonEnds) %>%
    mutate(
      exonStarts = strsplit(as.character(exonStarts), ","),
      exonEnds = strsplit(as.character(exonEnds), ",")
    ) %>%
    unnest(c(exonStarts, exonEnds)) %>%
    mutate(
      exonStarts = as.numeric(exonStarts),
      exonEnds = as.numeric(exonEnds)
    ) %>%
    group_by(row_id) %>%
    mutate(
      exon_idx = row_number()
    ) %>%
    ungroup()
  
  # Split and unnest blockStarts and blockEnds into long format
  block_df <- tx_domain_df %>%
    select(row_id, blockStarts, blockEnds) %>%
    mutate(
      blockStarts = strsplit(as.character(blockStarts), ","),
      blockEnds = strsplit(as.character(blockEnds), ",")
    ) %>%
    unnest(c(blockStarts, blockEnds)) %>%
    mutate(
      blockStarts = as.numeric(blockStarts),
      blockEnds = as.numeric(blockEnds)
    ) %>%
    group_by(row_id) %>%
    mutate(
      block_idx = row_number(),
      total_blocks = n(),
      is_single_block = total_blocks == 1,
      is_first_block = block_idx == 1,
      is_last_block = block_idx == total_blocks
    ) %>%
    ungroup()
  
  # Map blocks to exons with the specified criteria
  block_exon_df <- block_df %>%
    inner_join(exon_df, by = "row_id", relationship = "many-to-many") %>%
    filter(
      # Ensure blocks are within exon boundaries
      blockStarts >= exonStarts & blockEnds <= exonEnds
    ) %>%
    filter(
      # For single-block domains
      (is_single_block & blockStarts >= exonStarts & blockEnds <= exonEnds) |
        # For multiple-block domains
        (!is_single_block & (
          # First block: blockEnds matches exonEnds
          (is_first_block & blockEnds == exonEnds) |
            # Last block: blockStarts matches exonStarts
            (is_last_block & blockStarts == exonStarts) |
            # Middle blocks: exact match on both starts and ends
            (!is_first_block & !is_last_block & blockStarts == exonStarts & blockEnds == exonEnds)
        ))
    ) %>%
    group_by(row_id, block_idx) %>%
    summarise(
      exon_idx = list(unique(exon_idx)),
      .groups = 'drop'
    ) %>%
    mutate(
      exon_count = lengths(exon_idx)
    )
  
  # Identify domains with invalid blocks (blocks that do not map uniquely to one exon)
  invalid_blocks <- block_exon_df %>%
    filter(exon_count != 1)
  
  domains_with_invalid_blocks <- unique(invalid_blocks$row_id)
  
  # Valid blocks (blocks that map uniquely to one exon)
  valid_blocks <- block_exon_df %>%
    filter(exon_count == 1) %>%
    mutate(
      exon_idx = unlist(exon_idx)
    )
  
  # Collect exon indices per domain, ordered by block_idx
  exon_indices_per_domain <- valid_blocks %>%
    arrange(row_id, block_idx) %>%
    group_by(row_id) %>%
    summarise(
      exon_indices = list(exon_idx),
      .groups = 'drop'
    )
  
  # Total number of blocks per domain
  block_counts <- block_df %>%
    group_by(row_id) %>%
    summarise(
      total_blocks = n(),
      .groups = 'drop'
    )
  
  # Compile domain information
  domain_info <- block_counts %>%
    left_join(exon_indices_per_domain, by = "row_id") %>%
    mutate(
      valid_blocks = lengths(exon_indices),
      # Check for any invalid blocks or exon skipping
      has_issue = (row_id %in% domains_with_invalid_blocks) |
        (valid_blocks != total_blocks) |
        (
          # Check for exon skipping
          sapply(exon_indices, function(exons) {
            if (length(exons) <= 1) {
              return(FALSE)  # No exon skipping possible with one exon
            } else {
              return(!all(diff(exons) == 1))
            }
          })
        )
    )
  
  # Filter out domains with issues
  valid_domains <- domain_info %>%
    filter(!has_issue) %>%
    select(row_id)
  
  # Keep only the valid domains in the original dataframe
  tx_domain_df <- tx_domain_df %>%
    semi_join(valid_domains, by = "row_id") %>%
    select(-row_id)
  
  return(tx_domain_df)
}
