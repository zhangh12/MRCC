
find.mle <- function(par.tsls, rdata, edata, par.pos){

  sol <- optimx(par.tsls, nlogL, gr = nscore, method = 'ucminf',
                rdata = rdata, edata = edata, par.pos = par.pos)

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