#' Normalize data (median method)
#'
#' @description Normalize data by dividing by median and multiplying by mean of medians
#' @param data a data frame .
#' @param sample_names column identfiers or indices marking which columns to normalized
#'
#' @details This function peforms median-normalization on the requested columns of a data frame.

normalize_data <- function(data, sample_names) {
  if(is.matrix(data)) data <- as.data.frame(data)

  sample_medians <- apply(data[sample_names],2,median,na.rm = TRUE)
  mean_of_medians <- mean(sample_medians)

  normalized_data <- lapply(sample_names, function(sample_name) {
    data[[sample_name]] / sample_medians[sample_name] * mean_of_medians
  })
  data[sample_names] <- normalized_data

  data
  # colnames(normalized_data) <- paste0("normalized_",sample_names)
  # cbind(data,normalized_data)
}
