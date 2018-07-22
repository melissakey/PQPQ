library(PQPQ)

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

result <- pqpq(testdata2, sample_names = sample_names, data_type = "Manually annotated", manually_annotated_fields = column_ids)



###########################################################
#
# call function piece by piece
#
###########################################################

library(tidyverse)
library(magrittr)

processed_data <- preprocess_pqpq_input(testdata2, sample_names = sample_names, protein_subset = "gi|114657944")
column_ids <- processed_data$column_ids

data_list <- peptide_selection(processed_data$data,
  processed_data$column_ids)

output <- filter_peptides(data_list, column_ids, action = 'mark')

a <- "gi|114657944"
output %>%
  filter(Accessions == "gi|114657944") %>%
  select(Sequence, A, B, C)
