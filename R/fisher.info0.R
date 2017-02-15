
fisher.info0 <- function(par, rdata, edata, par.pos){

  id3 <- intersect(rownames(rdata$data), rownames(edata$data))
  id1 <- setdiff(rownames(rdata$data), id3)
  id2 <- setdiff(rownames(edata$data), id3)

  id10 <- setdiff(rownames(rdata$data[rdata$data[, rdata$vd] == 0, , drop = FALSE]), id3)
  id11 <- setdiff(rownames(rdata$data[rdata$data[, rdata$vd] == 1, , drop = FALSE]), id3)

  m1 <- length(id1)
  m2 <- length(id2)
  m3 <- length(id3)

  n10 <- length(id10)
  n11 <- length(id11)

  n0 <- n10 + m3
  n1 <- n11

  hess <- hessian0(par, rdata, edata, par.pos)
  J1 <- -hess[, 'a', drop = FALSE]
  info <- -hess - (n0 + n1)/n0 /n1 * (J1 %*% t(J1))
  info

}
