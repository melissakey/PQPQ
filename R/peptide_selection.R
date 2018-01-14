#' Select Peptides
#'
#' @description
#' This function does some additional processing, and validates the peptides assigned to each protein.
#' @param input_data A data frame or list of data frames
#' @param column_identifiers a list of labeled variable names.
#' @param the_limit "peptide sum intensity limit"
#' @param p_val "maximum p-value of the correlation coefficient
#' @param peptide_confidence_limit "high confidence limit", or the minimum peptide confidence score
#' @param do_normalization should the data be normalized?
#'
#' @details
#' If the data is entered as a list of data frames (1 per replicate), the data is split across replicates and labeled
#'
#' @return
#' A list, containing the cleaned data.frame and the list of labeled variable names.
#' @export




peptide_selection <- function(input_data, # data frame with input data
  column_identifiers,                     # a list, containing numeric or character vectors identifying the relevant columns of the data set.
  the_limit = 0,                          # "peptide sum intensity limit"
  p_val = 0.4,                            # maximum p-value of correlation coefficient
  peptide_confidence_limit = 95,          # "high confidence limit"
  # options
  do_normalization = TRUE,                # normalizes data to median
  # subsetting
  protein_subset = NULL                   # NULL, or list of proteins on which to perform filtering.
) {

  ######-------------------------------------------------------------------------
  #
  # Combine replicates into a single data frame.
  #
  ######-------------------------------------------------------------------------
  if (inherits(input_data,"list")) {
    replicates <- length(df_list)
    df <- plyr::ldply(seq_along(df_list), function(i) {
      within(df_list[[i]], {replicate <- LETTERS[i]})
    })
  } else if (inherits(input_data,"data.frame")) {
    replicates <- 1
    df <- input_data

  } else {
    stop("'input_data' must have class 'data.frame' or 'list'.")
  }





  ######-------------------------------------------------------------------------
  #
  # Data pre-processing
  #
  ######-------------------------------------------------------------------------

  # remove low-quality peptides:
  df[column_identifiers$sample_names] <- lapply(df[column_identifiers$sample_names], function(x, the_limit) {
    x[x < the_limit] <- NA
    x
  }, the_limit = the_limit)


  # normalization - if requested (separately for each replicate).
  if (do_normalization) {
    if(replicates > 1)
      df <- plyr::ddply(df, .(replicate), normalize_data, sample_names = column_identifiers$sample_names)
    else
      df <- normalize_data(df, column_identifiers$sample_names)

  }

  # replicates are already combined into a single data frame.


  ######-------------------------------------------------------------------------
  #
  # Peptide Selection
  #
  ######-------------------------------------------------------------------------

  # if (do_peptide_selection) {
  # peptide_selection_df <- df

  protein_list <- plyr::dlply(df, column_identifiers$protein_id, function(protein_df) {
    dl <- list(
      df = protein_df,
      matrix_data = t(as.matrix(protein_df[column_identifiers$sample_names])),
      confidence =protein_df[[column_identifiers$confidence]]
    )
    colnames(dl$matrix_data) <- names(dl$confidence) <- apply(protein_df[column_identifiers$peptide_ids],1,paste,collapse = "_")
    dl
  })

  result <- select_peptide_data(protein_list,
    data_matrix_name = "matrix_data",
    score_limit = score_limit,
    p_val =  p_val,
    peptide_confidence_limit = peptide_confidence_limit)
  # }
  class(result) <- "pqpq"
  result

}

