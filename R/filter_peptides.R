#' filter peptides
#'
#' Filter or mark non-correlating peptides.
#' @param pqpq_output a pqpq object
#' @param column_ids a list of labeled variable names.
#' @param action either "filter" or "mark".  See Details.
#' @details
#' This function takes an object of class pqpq and returns a data frame.  If \code{action} = 'filter',
#' peptides which do not meet the PQPQ criteria are discarded.  If \code{action} = 'mark', then a new column
#' \code{keep} is added which is marked TRUE for peptides retained by the algoithm.
#' @export
#'
#'


filter_peptides <- function(pqpq_obj,column_ids,action = c("filter","mark")) {

  # bookkeeping
  action <- match.arg(action)
  if(!inherits(pqpq_obj, 'pqpq')) stop("'pqpqp_output' must be of class 'pqpq'")

  plyr::ldply(pqpq_obj,
    function(x, column_ids) {
      n_clusters <- length(x$correlating_peptides)
      df <- x$df
      peptide_id <- apply(
        df[column_ids$peptide_ids], 1, paste, collapse = "_")
      if (n_clusters == 0) {
        if (action == 'filter'){
          df$keep <- FALSE
        } else {
          df <- NULL
        }
      } else {
        FILTER <- peptide_id %in% unique(unlist(x$correlating_peptides))
        if (action == 'filter') {
          df$keep <- FILTER
        } else df <- df[FILTER, ]
        names(x$correlating_peptides) <- LETTERS[1:n_clusters]
        members <- lapply(x$correlating_peptides,
          function(cp, peptide_id) {
            peptide_id %in% cp
          },
          peptide_id = peptide_id
        )
        df <- cbind(df,as.data.frame(members))
      }
      df
    },
    column_ids = column_ids
  )
}
