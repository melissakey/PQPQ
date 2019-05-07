#' select_peptide_data
#'
#' Select correlating peptides and subset them into clusters.
#' @param X the input data - see details.
#' @param data_matrix_name the name of the data matrix with the (normalized) peak intensities
#' @param score_limit currently unused
#' @param p_val the p-value associated with the significance of the correlation between 2 peptides
#' @param cut_point the inconsistency cut-off for clustering peptides by correlation distance.
#'
#'

select_peptide_data <- function(X, data_matrix_name, score_limit, p_val, peptide_confidence_limit,
  cut_point = 1) {

  ################################################
  #  R implementation of PQPQ
  #  by Jenny Forshed 2011-05-16
  #  port to R by Melissa Key
  ################################################

  result <- lapply(X, function(protein) {
    if(!(data_matrix_name %in% names(protein))) stop("'data_matrix_name' does not point to a valid matrix")

    confidence <- protein$confidence

    Nk <- ncol(protein[[data_matrix_name]])

    # create empty list to put warningings in
    protein_warnings <- vector("list", 5)

    # Do this once here, then pass the relevant results instead of recalculating it a
    # million times calculate correlation and associated p-values (p-value
    # calculation requires Hmisc package)

    r_cor <- Hmisc::rcorr(protein[[data_matrix_name]])


    # proteins supported by > 300 peptides select the 300 most confident peptides for
    # the actual protein
    if (Nk > 300) {
      index <- order(confidence,decreasing = TRUE)[1:300]
      protein_warnings[[4]] <- "protein supported by > 300 peptides"
    } else {
      index <- 1:Nk
    }
    if(Nk == 1) return(NULL)


    ##### select correlating peptides #####
    correlation_list <- lapply(seq(1, length(index) - 1), function(k) {
      index_pairs <- data.frame(k = index[k], l = index[-(1:k)])  # skipping protein/peptide IDs

      within(index_pairs, {
        p <- r_cor$P[k[1], l]
        r2 <- r_cor$r[k[1], l]
        conf_l <- confidence[l]
        conf_k <- confidence[k]
        conf <- (conf_k >= peptide_confidence_limit) & (conf_l >= peptide_confidence_limit)
        pep_l <- colnames(protein[[data_matrix_name]])[l]
        pep_k <- colnames(protein[[data_matrix_name]])[k]
      })
    })
    correlation_data <- do.call("rbind", correlation_list)
    # correlation_data <- plyr::ldply(seq(1, length(index) - 1), function(k) {
    #   index_pairs <- data.frame(k = index[k], l = index[-(1:k)])  # skipping protein/peptide IDs
    #
    #   within(index_pairs, {
    #     p <- r_cor$P[k[1], l]
    #     r2 <- r_cor$r[k[1], l]
    #     conf_l <- confidence[l]
    #     conf_k <- confidence[k]
    #     conf <- (conf_k >= peptide_confidence_limit) & (conf_l >= peptide_confidence_limit)
    #     pep_l <- colnames(protein[[data_matrix_name]])[l]
    #     pep_k <- colnames(protein[[data_matrix_name]])[k]
    #   })
    # })

    # Now we use subsetting of the above data to pick which pairs we actually want to use
    # first try: p-value < p_val and r2 > 0 and both confidences > peptide_confidence_limit
    index_selecting_proteins <- subset(correlation_data, p < p_val & r2 > 0 &
        conf, c(k, l))
    # if that didn't pick any up, relax the criteria on the p_value (not other criteria)
    if (nrow(index_selecting_proteins) == 0) {
      index_selecting_proteins <- subset(correlation_data, r2 > 0 & conf, c(k,
        l))
    }
    # if that doesn't work, throw up our hands and move on to the next protein.
    skip_rest <- FALSE
    if (nrow(index_selecting_proteins) == 0) {
      protein_warnings[[2]] <- "peptide conf too low to support protein ID or peptides don't correlate"
      skip_rest <- TRUE
    }

    if (!skip_rest){

      best_peptides <- unique(unlist(index_selecting_proteins))

      ##### Hierarchical clustering #####
      if (length(best_peptides) > 3) {
        # create dendrogram
        hclust_obj <- hclust(d = as.dist(1 - r_cor$r[best_peptides, best_peptides]), method = "single")

        # cut dendorgram into clusters
        class_assignment <- inconsistent(hclust_obj, cut_point = cut_point)$clusters
      } else {
        # otherwise, just put everything in one class.
        class_assignment <- data.frame(label = colnames(protein[[data_matrix_name]])[best_peptides], cluster = 1)
      }

      # sort classes by size.  this actually is important later!
      num_class_members <- with(class_assignment, tapply(cluster, cluster, length))



      # select the model peptide with the highest intensity (from each class)
      model_peptides_list <- split(class_assignment, class_assignment$cluster)
      model_peptides <- sapply(model_peptides_list, function(class, data) {
        names(which.max(colMeans(data[, as.character(class$label), drop = FALSE])))
      },data = protein[[data_matrix_name]])

      # model_peptides <- plyr::daply(class_assignment, plyr::.(cluster), function(class, data) {
      #   names(which.max(colMeans(data[, as.character(class$label), drop = FALSE])))
      # },data = protein[[data_matrix_name]])

      indices_for_peptides_correlating_with_model_peptides <- plyr::llply(model_peptides,
        function(peptide) {
          subset(correlation_data, (pep_k == peptide | pep_l == peptide) & p <
              p_val & r2 > 0)
        })
      correlating_peptides <- lapply(indices_for_peptides_correlating_with_model_peptides,
        with, {
          unique(c(pep_k, pep_l))
        })

      # order by the number of high-confidence peptides initially added to the cluster
      correlating_peptides <- correlating_peptides[order(num_class_members, decreasing = TRUE)]

      # loop through each pair of clusters.  If the clusters share over 50% of peptides
      # using the # of peptides in the 1st cluster as a reference, then remove the
      # second cluster.
      if (length(correlating_peptides) > 1) {
        end_clusters <- lapply(seq(2, length(correlating_peptides)), function(k) {
          overlap_degree <- sapply(seq(1, k - 1), function(l) {
            n_joint_peptides <- length(intersect(correlating_peptides[[k]], correlating_peptides[[l]]))
            n_peptides_in_k <- length(correlating_peptides[[k]])
            n_joint_peptides/n_peptides_in_k * 100
          })

          if (mean(overlap_degree) > 50)
            NULL else k
        })
        end_clusters <- c(1, unlist(end_clusters))
      } else {
        end_clusters <- 1
      }
      # create warning if more than one cluster is in the output.
      if(length(end_clusters) > 1) protein_warnings[[1]] <- 'non-overlapping peptide clusters'

      # create separate warning if a peptide is assigned to multiple clusters.
      # this warning is combined with warning #1 in original code.
      peptide_cluster_count <- table(unlist(correlating_peptides[end_clusters]))
      if(max(peptide_cluster_count) > 1) protein_warnings[[5]] <- 'peptides assigned to more than one cluster (overlapping clusters)'

      # attach the results to the original data contents.
      result <- c(protein,
        list(
          correlating_peptides = correlating_peptides[end_clusters],
          model_peptides = model_peptides,
          warnings = protein_warnings
        )
      )
    }
    else{
      result <- c(protein,
        list(
          warnings = protein_warnings
        )
      )
    }
    # class(result) <- 'pqpq'
    result
  })



}
