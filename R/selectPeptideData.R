selectPeptideData <- function(X,scoreLimit,pVal,peptideConfidenceLimit,correlationPvalue,dend.height = 1){
  require(plyr)
  require(Hmisc)

  ##################################
  # R implementation of PQPQ       #
  # by Jenny Forshed 2011-05-16    #
  # conversion by Melissa Key      #
  ##################################

  result <- ldply(X,with,{
    Nk <- ncol(data)

    # Do this once here, then pass the relevant results instead of recalculating it a million times
    # calculate correlation and associated p-values (p-value calculation requires Hmisc package)
    R.cor <- rcorr(data)


    # proteins supported by > 300 peptides
    # select the 300 most confident peptides for the actual protein
    # if(Nk > 300) index <- order(conf,decreasing = TRUE)[1:300]
    # else
    index <- 1:Nk
    protein.warnings <- vector("list",4)

    ##### select correlating peptides #####
    {
    # switch these if we cut off all but 300 peptides
    # correlationData <- ldply(2:length(index) -1,function(k){
    correlationData <- ldply(head(index,-1),function(k){

      within(data.frame(
        k = index[k],
        l = index[-(1:k)]
        # skipping protein/peptide IDs
      ),{
        p <- R.cor$P[k[1],l]
        r2 <- R.cor$r[k[1],l]
        Conf.l <- Conf[l]
        Conf.k <- Conf[k]
        conf <- (Conf.k >= peptideConfidenceLimit) & (Conf.l >= peptideConfidenceLimit)
        pep.l <- colnames(data[,l,drop = FALSE])
        pep.k <- colnames(data[,k,drop = FALSE])
      })
    })
    indexSelectedProteins = subset(correlationData,p < pVal & r2 > 0 & conf, c(k,l))
    if(is.null(indexSelectedProteins))
      indexSelectedProteins = subset(correlationData,r2 > 0 & conf, c(k,l))
    if()


    ##### Hierarchical clustering #####
    if(Nk > 3){
      all.ind <- unique(unlist(indexSelectedProteins))
      dend.obj <- as.dendrogram(hclust(d = as.dist(1 - R.cor$r[all.ind,all.ind]),method = "single"))
      clusters <- cut(dend.obj,h = dend.height)$lower
      # unique.classes <- length(clusters)
      class.assignment <- ldply(1:length(clusters),function(i){
        data.frame(
          label = labels(clusters[[i]]),
          cluster = i
        )
      })
    } else{
      class.assignment <- data.frame(
        label = colnames(data),
        cluster = 1
      )
    }
    noClassMembers <- with(class.assignment,tapply(cluster,cluster,length))


    # select the model peptide with the highest intensity (from each class)
    model.peptides <- daply(class.assignment,.(cluster),with,{
      names(which.max(colMeans(data[,as.character(label),drop = FALSE])))
    })

    indicesForPeptidesCorrelatingWithModelPeptides <- llply(model.peptides,function(peptide){
      subset(correlationData,
        (pep.k == peptide | pep.l == peptide) & p < pVal & r2 > 0)
    })
    correlating.peptides <- llply(indicesForPeptidesCorrelatingWithModelPeptides,with,{
      unique(c(pep.k,pep.l))
    })

  overlap.degree <- unlist(c(1,llply(seq(2, length(correlating.peptides)) ,function(k){
      overlap.degree <- laply(seq(1,k-1),function(l){
        n.joint.peps <- length(intersect(correlating.peptides[[k]],correlating.peptides[[l]]))
        n.k.peps <- length(correlating.peptides[[k]])
        n.joint.peps/n.k.peps * 100
      })
      if(overlap.degree > 50) NULL
      else k
    })))


  overlap.mat <- matrix(0,nrow = length(correlating.peptides),ncol = length(correlating.peptides))
  overlap.mat[lower.tri(overlap.mat)] <- unlist(overlap.degree)
  overlap.mat <- overlap.mat + t(overlap.mat)
    tmp <- unlist(overlap.degree)
    class(tmp) <- "dist"
    tmp2 <- as.matrix(tmp)

    result <- list(
      warnings = rep(0,4),

    )

  })



}
