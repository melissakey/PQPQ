#' PQPQ
#'
#' @description This function applies the PQPQ algorithm to proteomic data to filter out poorly correlating peptides and cluster the remaining peptides into likely proteoforms.
#' @param df The data frame with proteomic data.
#' @param sample_names A character vector identifying the columns holding the quantitative peptide data.  Optional IF (1) the data-type is not "manually annotated" and all columns are desired.
#' @param protein_subset The subset of proteins on which to apply the filter.
#' @param data_type One of "Protein Pilot", "Spectrum Mill", "Proteome Discoverer", or "Manually annotated" - Used to determine the columns containing different data elements.
#' @param normalize_data Should the data be normalized?  Default is TRUE for compatibility with the Matlab version.
#' @param correlation_p_value Correlations with p-values below this threshold are determined to be significant.
#' @param high_confidence_limit Minimum confidence for a peptide to be considered highly likley to be identified correctly.
#' @param peptide_sum_intensity_limit Minimum intensity for peptide to be included.
#' @param separate_multiple_protein_IDs If TRUE, peptides assigned to multiple proteins are copied and analyzed multiple times as part of each protein to which it is assigned.  If FALSE, peptides are assigned to the protein group, and analyzed once.  Default is FALSE.
#' @param manually_annotated_fields If data_type is "Manually annotated", this is a list of the columns needed to complete the filtering: protein_id, confidence, and peptide_ids.  See examples
#' @details
#' This function performs the PQPQ process - including preprocessing, peptide_selection, and filtering.
#' @export
#' @return A data frame identifying which peptides are kept, and which proteoforms they are assigned to.
#' @examples
#' data("testdata2")
#'
#' data("testdata2")
#' sample_names <- stringr::str_subset(names(testdata2), "^Area")[-7]
#' column_ids <- list(
#'   protein_id = "Accessions",
#'   confidence = 'Conf',
#'   peptide_ids = "Sequence"
#' )
#'
#' result <- pqpq(testdata2, sample_names = sample_names, data_type = "Manually annotated", manually_annotated_fields = column_ids)


pqpq <- function(data,                    # data frame with input data
  sample_names = NULL,
  protein_subset = NULL,
  data_type = c("Protein Pilot", "Spectrum Mill", "Proteome Discoverer", "Manually annotated"),
  normalize_data = TRUE,
  correlation_p_value = 0.4,
  high_confidence_limit = 95,
  peptide_sum_intensity_limit = 0,
  separate_multiple_protein_IDs = FALSE,
  manually_annotated_fields = NULL,
  action = c('mark', 'filter')
) {

  action <- match.arg(action)
  data_type <- match.arg(data_type)

  if(inherits(data, "data.frame")){
    processed_data <- preprocess_pqpq_input(data,
      data_type = data_type,
      protein_subset = protein_subset,
      separate_multiple_protein_IDs = separate_multiple_protein_IDs,
      sample_names = sample_names,
      manually_annotated_fields
    )
  } else if(inherits(data, "list")) {
    processed_data_list <- lapply(data, function(df) {
      preprocess_pqpq_input(df,
        data_type = data_type,
        protein_subset = protein_subset,
        separate_multiple_protein_IDs = separate_multiple_protein_IDs,
        sample_names = sample_names,
        manually_annotated_fields
      )
    })

    sample_names_tmp <- lapply(processed_data_list, function(x) x$column_ids$sample_names)
    if(length(unique(sample_names_tmp)) == 1) {
      column_ids <- processed_data_list[[1]]$column_ids
    } else {
      column_ids <- unique(unlist(sample_names_tmp))
    }


    processed_data <- list(
      data = lapply(processed_data_list, with, data),
      column_ids = column_ids
    )
  }

  data_list <- with(processed_data, peptide_selection(data,
    column_ids, the_limit = peptide_sum_intensity_limit, p_val = correlation_p_value, peptide_confidence_limit = high_confidence_limit, do_normalization = normalize_data))

  pqpq_warnings <- ldply(data_list, function(x) as.data.frame(x$warnings))
  output <- filter_peptides(data_list, processed_data$column_ids, action = action)
  output <- list(
    result = output,
    warnings = pqpq_warnings
  )
  class(output) <- 'pqpq'
}

