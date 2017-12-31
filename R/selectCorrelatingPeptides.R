selectCorrelatingPeptides <- function(data, Conf, index, peptideConfidenceLimit,correlationPvalue){
  require(Hmisc)
  ##################################
  # R implementation of PQPQ       #
  # by Jenny Forshed 2011-05-16    #
  # conversion by Melissa Key      #
  ##################################

  # delete rows with zero intensities
  # commented out in original code
  # min.row.val <- apply(data,1,min)

  # find a peptide pattern by checking if the data from different
  # peptides are correlating, include only peptides with peptide confidence >
  # peptideConfidenceLimit

  # ORDER CHANGE FROM ORIGINAL CODE
  # first figure out which correlations we actually want
  indexCL <- index[Conf[index] >= peptideConfidenceLimit]

  # calculate correlations only for that subset
  R.cor <- rcorr(data[,indexCL])

  correlationData <- ldply(2:length(indexCL) -1,function(k){
    within(data.frame(
      k = indexCL[k],
      l = indexCL[-(1:k)]
      # skipping protein/peptide IDs
    ),{
      p <- R.cor$P[k[1],l]
      r2 <- R.cor$r[k[1],l]
    })
  })
  list(
    data = data,
    correlationData = correlationData,

  )
}
