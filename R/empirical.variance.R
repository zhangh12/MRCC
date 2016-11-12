
empirical.variance <- function(par, rdata, edata, par.pos){

  h <- hessian(par, rdata, edata, par.pos)
  i <- fisher.info(par, rdata, edata, par.pos)

  solve(h) %*% i %*% solve(h)

}
