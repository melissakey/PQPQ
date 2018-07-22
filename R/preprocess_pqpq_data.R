#' Preprocess data for entry into PQPQ software
#'
#' @description
#' This function is used to prepare Protein Pilot data for analysis in PQPQ.
#' @param df A data frame (see details)
#' @param data_type one of 'Protein Pilot', 'Spectrum Mill', or 'Proteome Discoverer', depending on the file origin.
#' @param protein_subset identifiers for the proteins to be analyzed.
#' @param separate_multiple_protein_IDs If true, peptides assigned to multiple proteins are duplicated for each protein to which they are assigned.  Otherwise, they are analyzed as a group.
#' @param sample_names A character vector identifying the columns holding sample data.
#' @param manually_annotated_fields A list containing annotation data.
#'
#' @details
#' Note: this is only tested on the test output from Protein Pilot provided by Jenny Forshed in the PQPQ package, and entered into R using the \code{link{openxlsx}} package.  The formatting of data from other
#'
#' @export
#' @return
#' A list, containing the cleaned data.frame and the list of labeled variable names.
#' @examples
#' # the column identifiers column - minimal example
#' data("testdata2")
#'
#' # minmal column_ids example
#' column_ids <- list(
#'    protein_id = "Accessions",
#'    confidence = "Conf",
#'    samples_names = grep("Area",names(testdata2),value = TRUE),
#'    peptide_ids = c("Sequence", "peptide_id")
#' )


preprocess_pqpq_input <- function(df,
  data_type = "Protein Pilot",
  protein_subset = NULL,
  separate_multiple_protein_IDs = FALSE,
  sample_names = NULL,
  manually_annotated_fields = NULL
) {

  df <- within(df,{
    if(!exists("peptide_id",inherits = FALSE)) {
      peptide_id <- 1:nrow(df)
    } else if(length(unique(peptide_id)) != length(peptide_id)) {
      peptide_id <- paste(peptide_id, 1:nrow(df), sep = "_")
    }
  })

  column_ids <- switch(data_type,
    `Protein Pilot` = {
      list(protein_id = "Accessions",
        protein_name = "Names",
        confidence = 'Conf',
        sample_names = {if(!is.null(sample_names)) sample_names else grep("Area",names(df), value = TRUE)},
        peptide_ids = c("Sequence","peptide_id")
      )
    },
    `Spectrum Mill` = {
      warning("check these results - Melissa was not able to test this code.")
      list(
        protein_id = 'accessionnumber',
        protein_name = 'entryname',
        confidence = "score",
        sample_names = {if(!is.null(sample_names)) sample_names else grep("iTRAQ_1",names(df),value = TRUE)},
        peptide_ids = c("sequence","peptide_id")
      )
    },
    `Proteome Discoverer` = {
      warning("check these results - Melissa was not able to test this code.")
      list(
        protein_id = c('ProteinAccessions'),
        confidence = {if(exists("df$IonScore")) "IonScore" else "XCorr"},
        sample_names = {if(!is.null(sample_names)) sample_names else grep("^11X",names(df),value = TRUE)},
        peptide_ids = c("sequence","peptide_id")
      )
    },
    `Manually annotated` = {
      if(is.null(manually_annotated_fields)) stop("Manual annotation requires `manually_annotated_fields` to be defined")

      manually_annotated_fields$sample_names = sample_names
      if(
        !is.list(manually_annotated_fields) |
          length(setdiff(c("protein_id", "confidence", "sample_names", "peptide_ids"), names(manually_annotated_fields) )) > 0)
        stop("Manual annotation requires `manually_annotated fields to be a list with elements `protein_id`, `confidence`, and `peptide_ids`; and `sample_names` must be defined.")

      manually_annotated_fields <- within(manually_annotated_fields, {
        peptide_ids <- unique(peptide_ids, "peptide_id")
      })
    },
    stop("Please give a valid data_type or manually list column names.")
  )
  # check data parameters
  column_ids <- check_column_names(df, column_ids)

  # filter out peptides
  if (data_type == "Protein Pilot") {
    df <- df[df$Annotation != 'auto - shared MS/MS', ]
  } else if (data_type == 'Proteome Discoverer') {
    warning("check to make sure correct proteins are used.")
    df <- subset(df, QuanUsage == 'Used' & ProteinGroups != 0)
  }

  if(separate_multiple_protein_IDs) {
    df <- plyr::ddply(df, column_ids$protein_id, function(protein_df) {
      protein_ids <- strsplit(protein_df[[column_ids$protein_id]][1], split = "; ")[[1]]

      rep_data <- plyr::ldply(protein_ids,
        function(split_ids, local_df) {
          local_df$protein_id <- split_ids
          local_df
        },
        local_df = protein_df
      )
    })
    column_ids <- within(column_ids,{
      protein_original_id <- protein_id
      protein_id <- "protein_id"
    })
  }

  if(!is.null(protein_subset)){
    list(
      column_ids = column_ids,
      data = df[df[[column_ids$protein_id]] %in% protein_subset, ]
    )

  } else {
    list(
      column_ids = column_ids,
      data = df
    )
  }


}
