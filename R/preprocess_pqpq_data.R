#' Preprocess data for entry into PQPQ software
#'
#' @description
#' This function is used to prepare Protein Pilot data for analysis in PQPQ.
#' @param df A data frame (see details)
#' @param data_type one of 'Protein Pilot', 'Spectrum Mill', or 'Proteome Discoverer', depending on the file origin.
#' @param protein_subset identifiers for the proteins to be analyzed.
#' @param separate_multiple_protein_IDs If true, peptides assigned to multiple proteins are duplicated for each protein to which they are assigned.  Otherwise, they are analyzed as a group.
#' @param column_ids A list of variables identifiers.
#'
#' @details
#' This function handles most of the pre-processing done in the PQPQ function including:
#' \itemize{
#' \item{If \code{column_ids} is NULL, creates a list of which column names identifying the relevant portions of the data set (confidence, protein ID, peptide ID, and samples).
#' This is subject to error depending on how the data was input into R and the data type.}
#' \item{confirms that the specified columns are in the data set.}
#' \item{removes peptides which do not meet PQPQ criterion. (in Protein Pilot, this includes proteins with annotation 'auto - shared MS/MS')}
#' \item{if separate_multiple_proteins_IDs = TRUE, assigns peptides belonging to multiple proteins to ALL listed proteins separately, not the protein group.}
#' \item{if a protein subset is given, filters out all other proteins.}
#'}
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
  column_ids = NULL
) {

  df <- within(df,{
    peptide_id <- 1:nrow(df)
  })

  if (is.null(column_ids)) {
    column_ids <- switch(data_type,
      `Protein Pilot` = {
        list(protein_id = "Accessions",
          protein_name = "Names",
          confidence = 'Conf',
          sample_names = grep("Area",names(df),value = TRUE),
          peptide_ids = c("Sequence","peptide_id")
        )
      },
      `Spectrum Mill` = {
        warning("check these results - Melissa was not able to test this code.")
        list(
          protein_id = 'accessionnumber',
          protein_name = 'entryname',
          confidence = "score",
          sample_names = grep("iTRAQ_1",names(df),value = TRUE),
          peptide_ids = c("sequence","peptide_id")
        )
      },
      `Proteome Discoverer` = {
        warning("check these results - Melissa was not able to test this code.")
        list(
          protein_id = c('ProteinAccessions'),
          confidence = {if(exists("df$IonScore")) "IonScore" else "XCorr"},
          sample_names = grep("^11X",names(df),value = TRUE),
          peptide_ids = c("sequence","peptide_id")
        )
      },
      stop("Please give a valid data_type or manually list column names.")
    )
  }
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
      parameters = column_ids,
      data = df[df[[column_ids$protein_id]] %in% protein_subset, ]
    )

  } else {
    list(
      column_ids = column_ids,
      data = df
    )
  }


}
