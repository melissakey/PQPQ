library(PQPQ)
library(tidyverse)
library(magrittr)
filter <- dplyr::filter


###########################################################
#
# best way to call the function
#
###########################################################


data("testdata2")
sample_names <- grep("^Area", names(testdata2), value = TRUE)[-7]
column_ids <- list(
  protein_id = "Accessions",
  confidence = 'Conf',
  peptide_ids = "Sequence"
)

output <- pqpq(testdata2, sample_names = sample_names)
output %>%
  filter(Accessions == "gi|114657944") %>%
  select(Sequence, A, B, C)

###########################################################
#
# test of multiple replicates
#
###########################################################


data("testdata2")
sample_names <- grep("^Area", names(testdata2), value = TRUE)[-7]
column_ids <- list(
  protein_id = "Accessions",
  confidence = 'Conf',
  peptide_ids = "Sequence"
)

output <- pqpq(list(A = testdata2, B = testdata2, C = testdata2), sample_names = sample_names)
output %>%
  filter(Accessions == "gi|114657944") %>%
  select(Sequence, A, B, C)




###########################################################
#
# call function piece by piece
#
###########################################################


processed_data <- preprocess_pqpq_input(testdata2, sample_names = sample_names)
column_ids <- processed_data$column_ids

data_list <- peptide_selection(processed_data$data,
  processed_data$column_ids)

output <- filter_peptides(data_list, column_ids, action = 'mark')
output %>%
  filter(Accessions == "gi|114657944") %>%
  select(Sequence, A, B, C)
