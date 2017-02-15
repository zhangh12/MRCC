
# given bet, find MLE of the rest of parameters
find.LRT.mle <- function(rdata, edata, bet){

  tsls <- find.LRT.tsls(rdata, edata, bet)

  ap <- align.parameter(tsls, null = TRUE)
  par.tsls <- ap$par
  par.pos <- ap$par.pos

  sol <- optimx(par.tsls, nlogL.LRT, gr = nscore.LRT, method = 'ucminf',
                rdata = rdata, edata = edata, par.pos = par.pos, bet = bet)

  # sol <- optimx(par.tsls, nlogL.LRT, gr = nscore.LRT, hess = nhessian.LRT, method = 'ucminf',
  #               rdata = rdata, edata = edata, par.pos = par.pos, bet = bet)

  if(!sol['convcode']){
    par <- as.vector(t(sol[1, names(par.tsls)]))
    names(par) <- names(par.tsls)

    if(check.LRT.mle(par, rdata, edata, par.pos, bet)){
      return(par)
    }
  }

  stop('optimx fails')

}
