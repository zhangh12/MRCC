
fisher.info0 <- function(par, rdata, edata, omega, par.pos){

  id10 <- omega$id10
  id11 <- omega$id11
  id30 <- omega$id30
  id31 <- omega$id31

  m10 <- length(id10)
  m11 <- length(id11)
  m30 <- length(id30)
  m31 <- length(id31)

  n0 <- m10 + m30
  n1 <- m11 + m31

  h0 <- hessian0(par, rdata, edata, omega, par.pos)
  J1 <- -h0[, 'a', drop = FALSE]
  info <- -h0 - (n0 + n1)/n0 /n1 * (J1 %*% t(J1))
  info

}
