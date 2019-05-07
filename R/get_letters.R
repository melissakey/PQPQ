make.chars <- function(length.out, case, n.char = NULL) {
  if(is.null(n.char))
    n.char <- max(1, ceiling(log(length.out, 26)))
  m <- sapply(n.char:1, function(x) {
    rep(rep(1:26, each = 26^(x-1)) , length.out = length.out)
  })
    m.char <- switch(case,
      'lower' = letters[m],
      'upper' = LETTERS[m]
    )
    m.char <- LETTERS[m]
    if(is.null(dim(m))) dim(m.char) <- c(1,n.char)
    else dim(m.char) <- dim(m)

    apply(m.char, 1, function(x) paste(x, collapse = ""))
}


get.letters <- function(length.out, case = 'upper'){
  max.char <- max(1, ceiling(log(length.out, 26)))
  grp <- rep(1:max.char, 26^(1:max.char))[1:length.out]
  unlist(lapply(unique(grp), function(n) make.chars(length(grp[grp == n]), case = case, n.char = n)))
}
