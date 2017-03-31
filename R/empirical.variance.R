
# generally over-estimate the variance
empirical.variance <- function(par, rdata, edata, omega, par.pos){

  h <- hessian(par, rdata, edata, par.pos)
  i <- fisher.info(par, rdata, edata, omega, par.pos)

  sqrt(diag(solve(h) %*% i %*% solve(h)))

}
