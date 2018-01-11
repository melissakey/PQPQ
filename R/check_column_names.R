#' Check column names
#'
#' @description
#' This function checks to make sure the variables required by PQPQ are in the data set.
#' @param df A data frame (see details)
#' @param column_identifiers A list of parameters in the data set.  See details.
#'
#' @details
#' This is a utility function designed to make sure that each parameter included in the column identfiers is present in the data set.
#' The parameter \code{column_identifer} should be a list with fields \code{protein_id} identifying a single varible with the proteins/Accessions, \code{confidence} identifying a single variable containing the  confidence for each peptide identity, \code{sample_names} identifying which columns contain the quantification data for each biological sample, and \code{peptide_ids} containing the columns uniquely identifying each peptide.
#'
#' @return
#' This function returns the parameter list, with any factors replaced by characters.
#' The fuction will currently stop if an error is detected, so that the parameters causing problems can be fixed by the user.
#'
#' @export



check_column_names <- function(df, column_identifiers) {

  ######-------------------------------------------------------------------------
  #
  # Check column identifiers
  #
  ######-------------------------------------------------------------------------

  # convert any factors to characters:
  column_identifiers[sapply(column_identifiers, is.factor)] <- lapply(
    column_identifiers[sapply(column_identifiers, is.factor)],
    as.character
  )

  # check to make sure required input variables are included.
  names_to_check <- c("confidence", "protein_id", "sample_names", "peptide_ids")
  for (var_name in names_to_check)
    if (!(var_name %in% names(column_identifiers))) stop(paste("'column_identifiers' must include", var_name))

  # check length where relevant
  if(length(column_identifiers$confidence) != 1) stop("'confidence' parameter must be of length 1, identifying the column containing confidence values.")
  if(length(column_identifiers$protein_id) != 1) stop("'protein_id' parameter must be of length 1, identifying the column containing protein identifiers.")


  # check each for invalid , going down the list:
  plyr::l_ply(names(column_identifiers), function(input) {
    input_vector <- column_identifiers[[input]]

    # convert each input to numeric identifier.
    col_index <- tryCatch(
      col_index <- {
        if(is.numeric(input_vector)) match(names(df)[input_vector], names(df))
        else match(input_vector, names(df))
      },
      error = function(cond) {
        message(paste("Could not obtain", input, "from given indices."))
      },
      warning = function(cond) {
        message(paste("A problem occurred obtaining", input,"from indices.  Original warning:"))
        message(cond)
      }
    )
    if (any(is.na(col_index))) {
      stop(paste("Could not obtain", input, "from given indices."))
    }
  })
  column_identifiers
}
