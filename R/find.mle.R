
find.mle <- function(par.tsls, rdata, edata, par.pos){

  t1 <- try(sol <- optimx(par.tsls, nlogL, gr = nscore, method = 'ucminf', hessian = TRUE, control = list(dowarn = FALSE),
                          rdata = rdata, edata = edata, par.pos = par.pos))

  if(('try-error' %in% class(t1)) | sol['convcode']){
    lower <- rep(-Inf, length(par.tsls))
    names(lower) <- names(par.tsls)
    lower['c'] <- 1e-6
    sol <- optimx(par.tsls, nlogL, gr = nscore, lower = lower, method = 'L-BFGS-B', hessian = TRUE, control = list(dowarn = FALSE),
                  rdata = rdata, edata = edata, par.pos = par.pos)
  }

  # sol <- optimx(par.tsls, nlogL, gr = nscore, hess = nhessian, method = 'ucminf',
  #               rdata = rdata, edata = edata, par.pos = par.pos)

  if(!sol['convcode']){
    par <- as.vector(t(sol[1, names(par.tsls)]))
    names(par) <- names(par.tsls)

    if(check.mle(par, rdata, edata, par.pos)){
      return(par)
    }
  }

  stop('optimx fails')

}