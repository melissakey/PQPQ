
###  Helper function to search through hclust-output and returns all members of the selected cluster.
find_nodes <- function(left_side, right_side, index) {
  if(left_side[index] < 0) {
    left_val <- abs(left_side[index])
  } else {
    left_val <- find_nodes(left_side, right_side, left_side[index])
  }
  if(right_side[index] < 0) {
    right_val <- abs(right_side[index])
  } else {
    right_val <- find_nodes(left_side, right_side, right_side[index])
  }
  c(left_val, right_val)
}

# vectorized function to calculate variance for each row of a matrix
# from https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
rowVar <- function(x,na.rm = FALSE) {
  rowSums((x - rowMeans(x,na.rm = na.rm))^2,na.rm = na.rm)/(rowSums(!is.na(x)) - 1)
}


