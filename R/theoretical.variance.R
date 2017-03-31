
# generally under-estimate the variance
theoretical.variance <- function(par, rdata, edata, omega, par.pos){

  h0 <- hessian0(par, rdata, edata, omega, par.pos)
  i0 <- fisher.info0(par, rdata, edata, omega, par.pos)

  sqrt(diag(solve(h0) %*% i0 %*% solve(h0)))

}
