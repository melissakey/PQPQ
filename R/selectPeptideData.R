#' select_peptide_data
#'
#' Select correlating peptides and subset them into clusters.
#' @param X the input data - see details.
#' @param score_limit currently unused
#' @param p_val the p-value associated with the significance of the correlation between 2 peptides
#' @param dend_height the cut-off for clustering peptides by correlation distance.
#' @keywords proteomics
#'
#'

select_peptide_data <- function(X, score_limit, p_val, peptide_confidence_limit,
    dend_height = 1) {


    ################################################
    #  R implementation of PQPQ
    #  by Jenny Forshed 2011-05-16
    #  port to R by Melissa Key
    ################################################

      result <- plyr::ldply(X, with, {
        Nk <- ncol(data)

        # Do this once here, then pass the relevant results instead of recalculating it a
        # million times calculate correlation and associated p-values (p-value
        # calculation requires Hmisc package)
        r_cor <- Hmisc::rcorr(data)


        # proteins supported by > 300 peptides select the 300 most confident peptides for
        # the actual protein if (Nk > 300) index <- order(conf,decreasing = TRUE)[1:300]
        # else
        index <- 1:Nk
        protein_warnings <- vector("list", 4)

        ##### select correlating peptides #####

        # switch these if we cut off all but 300 peptides correlationData <-
        # ldply(2:length(index) -1,function(k){
        correlation_data <- plyr::ldply(head(index, -1), function(k) {
            index_pairs <- data.frame(k = index[k], l = index[-(1:k)]  # skipping protein/peptide IDs
)
            within(index_pairs, {
                p <- r_cor$P[k[1], l]
                r2 <- r_cor$r[k[1], l]
                conf_l <- Conf[l]
                conf_k <- Conf[k]
                conf <- (conf_k >= peptide_confidence_limit) & (conf_l >= peptide_confidence_limit)
                pep_l <- colnames(data[, l, drop = FALSE])
                pep_k <- colnames(data[, k, drop = FALSE])
            })
        })

        index_selecting_proteins <- subset(correlation_data, p < p_val & r2 > 0 &
            conf, c(k, l))

        if (is.null(index_selecting_proteins)) {
            index_selecting_proteins <- subset(correlation_data, r2 > 0 & conf, c(k,
                l))
        }
        if (is.null(index_selecting_proteins)) {
            protein_warnings[[2]] <- "peptide conf too low to support protein ID or peptides don't correlate"
        }


        ##### Hierarchical clustering #####
        if (Nk > 3) {
            all_ind <- unique(unlist(index_selecting_proteins))
            dend_obj <- as.dendrogram(hclust(d = as.dist(1 - r_cor$r[all_ind, all_ind]),
                method = "single"))
            clusters <- cut(dend_obj, h = dend_height)$lower
            # unique.classes <- length(clusters)
            class_assignment <- ldply(1:length(clusters), function(i) {
                data.frame(label = labels(clusters[[i]]), cluster = i)
            })
        } else {
            class_assignment <- data.frame(label = colnames(data), cluster = 1)
        }
        noClassMembers <- with(class_assignment, tapply(cluster, cluster, length))


        # select the model peptide with the highest intensity (from each class)
        model_peptides <- daply(class_assignment, .(cluster), with, {
            names(which.max(colMeans(data[, as.character(label), drop = FALSE])))
        })

        indices_for_peptides_correlating_with_model_peptides <- plyr::llply(model_peptides,
            function(peptide) {
                subset(correlationData, (pep_k == peptide | pep_l == peptide) & p <
                  p_val & r2 > 0)
            })
        correlating_peptides <- plyr::llply(indices_for_peptides_correlating_with_model_peptides,
            with, {
                unique(c(pep_k, pep_l))
            })
        correlating_peptides <- correlating_peptides[order(sapply(correlating_peptides,
            length, decreasing = TRUE))]

        # loop through each pair of clusters.  If the clusters share over 50% of peptides
        # using the # of peptides in the 1st cluster as a reference, then remove the
        # second cluster.
        overlap_degree <- plyr::llply(seq(2, length(correlating_peptides)), function(k) {
            overlap_degree <- plyr::laply(seq(1, k - 1), function(l) {
                n_joint_peptides <- length(intersect(correlating_peptides[[k]], correlating_peptides[[l]]))
                n_peptides_in_k <- length(correlating_peptides[[k]])
                n_joint_peptides/n_peptides_in_k * 100
            })
            if (overlap_degree > 50)
                NULL else k
        })
        overlap_degree <- c(1, unlist(overlap_degree))

        # overlap.mat <- matrix(0,nrow = length(correlating.peptides),ncol =
        # length(correlating.peptides)) overlap.mat[lower.tri(overlap.mat)] <-
        # unlist(overlap.degree) overlap.mat <- overlap.mat + t(overlap.mat) tmp <-
        # unlist(overlap.degree) class(tmp) <- 'dist' tmp2 <- as.matrix(tmp)

        # result <- list( warnings = rep(0,4), )

    })



}
