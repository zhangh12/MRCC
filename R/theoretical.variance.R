
theoretical.variance <- function(par, rdata, edata, par.pos){

  h0 <- hessian0(par, rdata, edata, par.pos)
  i0 <- fisher.info0(par, rdata, edata, par.pos)

  solve(h0) %*% i0 %*% solve(h0)

}
