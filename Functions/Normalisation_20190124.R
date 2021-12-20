# Thu Jan 24 13:21:55 2019 ------------------------------
# Normalisation functions pqn ---------------------------------------------


# x <- quant.data.3$data[[1]]
# 1. pqn KODAMA package ---------------------------------------------------
library(KODAMA)
pqn_kodama_new <- function(x) {
  y <- dcast(data = x,Sample~Name,value.var = "Value")
  z <- column_to_rownames(.data = y,var = "Sample")
  z[is.na(z)] <- 0
  a=z[ ,which(colSums(z) != 0)]
  ref <- apply(X = a,MARGIN = 2,FUN = function(x) median(x))
  b <- normalization(a,method = "pqn",ref = ref)$newXtrain
  c <- rownames_to_column(.data = b,var = "Sample")
  d <- melt(c,id.vars = "Sample")
  colnames(d)[2:3] <- c("Name","Value")
  as_tibble(d)
}


# 2. Selfwritten_pqn_function ---------------------------------------------

pqn_PL_MS1_data <- function(x) {
  y <-  pivot_wider(data=x, id_cols = c("Sample"),names_from = lipid,values_from= Value)
  z <- column_to_rownames(.data = y,var = "Sample")
  z[is.na(z)] <- 0
  a=z[ ,which(colSums(z) != 0)]
  ref <- apply(X = a,MARGIN = 2,FUN = function(x) median(x))
  b <- matrix(nrow=nrow(a),ncol=ncol(a))
  row.names(b) <- row.names(a)
  colnames(b) <- colnames(a)
  for (i in seq_len(nrow(a))) {
    b[i,] <- local({
      a <- unlist(a[i,])
      b <- c(unlist(ref))
      c <- a/b
      d <- median(c,na.rm = T)
      a/d})
  }
  c <- rownames_to_column(.data = as.data.frame(b),var = "Sample")
  d <- pivot_longer(data=c, id_cols = c("Sample"),names_to = "lipid",values_to = "Value")
  as_tibble(d)
}
