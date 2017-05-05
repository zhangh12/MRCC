
working.variance <- function(par, rdata, edata, omega, par.pos){

  h <- hessian(par, rdata, edata, par.pos)
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

  v <- solve(-h)
  v['a', 'a'] <- v['a', 'a'] - (n0+n1)/n0/n1
  v

}
