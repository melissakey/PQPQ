#' filter peptides
#'
#' Filter or mark non-correlating peptides.
#' @param pqpq_output a pqpq object
#' @param data_matrix_name the name of the data matrix with the (normalized) peak intensities
#' @param score_limit currently unused
#' @param p_val the p-value associated with the significance of the correlation between 2 peptides
#' @param dend_height the cut-off for clustering peptides by correlation distance.
#' @export
#'
#'


filter_peptides <- function(pqpq_output,column_ids) {
  if(!inherits(pqpqp_output))
  plyr::ldply(pqpq_output,
    function(x, column_ids) {
      n_clusters <- length(x$correlating_peptides)
      peptide_id <- apply(x$df[column_ids$peptide_ids], 1, paste,collapse = "_")
      if (n_clusters == 0) {
        df <- NULL
      } else {
        FILTER <- peptide_id %in% unique(unlist(x$correlating_peptides))
        df$keep <- FILTER
        if (n_clusters > 1) {
          names(x$correlating_peptides) <- LETTERS[1:n_clusters]
          members <- lapply(x$correlating_peptides,
            function(cp, peptide_id) {
              peptide_id %in% cp
            },
            peptide_id = peptide_id
          )
          df <- cbind(df,as.data.frame(members))
        }
      }
      df
    },
    column_ids = column_ids
  )
}
