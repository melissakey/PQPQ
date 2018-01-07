normalize_data <- function(data, sample_names) {

  sample_medians <- apply(data,2,median,na.rm = TRUE)
  mean_of_medians <- mean(sample_medians)

  normalized_data <- plyr::llply(sample_names, function(sample_name) {
    data[, sample_name] / sample_medians[sample_name] * mean_of_medians
  })
  data[, sample_names] <- normalized_data[, sample_names]

  # colnames(normalized_data) <- paste0("normalized_",sample_names)
  # cbind(data,normalized_data)
}
