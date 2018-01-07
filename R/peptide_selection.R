peptide_selection <- function(df_list,                    # data frame with input data
  input_parameters,                     # a list, containing numeric or character vectors identifying the relevant columns of the data set.
  # column identifiers
  # # sample_names,                         # numeric or character vector of columns with AUC data
  # # confidence,                           # numeric or character vector of column with peptide confidence level
  # # protein_ids,                           # numeric or character vector of column with protein identity ("Accessions")
  # # peptide_ids,                           # numeric or character vector of columns with peptide identifier (e.g. sequence/charge)
  # parameters
  the_limit = 0,                      # "peptide sum intensity limit"
  p_val = 0.4,                          # maximum p-value of correlation coefficient
  peptide_confidence_limit = 95,         # "high confidence limit"
  dend_height = 1,                          # height of dendrogram at which to apply split.
  # options
  do_normalization = TRUE,              # normalizes data to median
  # subsetting
  protein_subset = NULL                 # NULL, or list of proteins on which to perform filtering.
) {

  ######-------------------------------------------------------------------------
  #
  # Combine replicates into a single data frame.
  #
  ######-------------------------------------------------------------------------

  df <- ldply(seq_along(df_list), function(i) {
    within(df_list[[i]],{
      replicate <- LETTERS[i]
    })
  })

  ######-------------------------------------------------------------------------
  #
  # Check column identifiers -- remove if handled earlier
  #
  ######-------------------------------------------------------------------------

  # convert any factors to characters:
  input_parameters[sapply(input_parameters, is.factor)] <- lapply(
    input_parameters[sapply(input_parameters, is.factor)],
    as.character
  )

  # check to make sure required input variables are included.
  names_to_check <- c("confidence", "protein_ids", "sample_names", "peptide_ids")
  for (var_name in names_to_check)
    if (!(var_name %in% names(input_variable))) stop(paste("input_variable must include", var_name))

  # check length where relevant
  if(length(input_parameters$confidence) != 1) stop("'confidence' parameter must be of length 1, identifying the column containing confidence values.")
  if(length(input_parameters$protein_ids) != 1) stop("'protein_ids' parameter must be of length 1, identifying the column containing protein identifiers.")


  # check each for invalid , going down the list:
  col_index <- plyr::llply(names(input_parameters), function(input) {
    input_vector <- input_parameters[[input]]

    # convert each input to numeric identifier.
    col_index <- tryCatch(
      col_index <- {
        if(is.numeric(input_vector)) match(names(df)[input_vector], names(df))
        else match(input_vector, names(df))
      }
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


  ######-------------------------------------------------------------------------
  #
  # Data pre-processing
  #
  ######-------------------------------------------------------------------------

  # remove low-quality peptides:
  df[col_index$sample_names] <- lapply(df[col_index$sample_names], function(x, the_limit) {
    x[x < the_limit] <- NA
    x
  }, the_limit = the_limit)


  # normalization - if requested (separately for each replicate).
  if (do_normalization) {
    df <- ddply(df, .(replicate), normalize_data, sample_names = sample_names)
  }

  # replicates are combined into a single data frame using this method.



  ######-------------------------------------------------------------------------
  #
  # Peptide Selection
  #
  ######-------------------------------------------------------------------------

  if (do_peptide_selection) {
    peptide_selection_df <- df
    names(peptide_selection_df[col_index$protein_ids]) <- "protein_id"

    protein_list <- dlply(df, .(protein_id), function(protein_df) {
      split_data <- lapply(col_index, function(x) {
        if (x == 'protein_ids') protein_df$protein_id[1]
        else
          t(as.matrix(protein_df[x]))
      })
    })

    result <- select_peptide_data(protein_list, sample_names, score_limit, p_val, peptide_confidence_limit,
      dend_height = 1)
  }


}

