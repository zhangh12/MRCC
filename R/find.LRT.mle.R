
find.LRT.mle <- function(rdata, edata, par.pos, bet, ncut = 50){

  tsls <- find.LRT.tsls(rdata, edata, bet)

  ap <- align.parameter(tsls)
  par.tsls <- ap$par
  par.pos <- ap$par.pos

  t1 <- try(par <- newton.raphson.LRT(par.tsls, rdata, edata, par.pos), silent = TRUE)
  if(('try-error' %in% class(t1)) || !check.LRT.mle(par, rdata, edata, par.pos)){
    s0 <- score.LRT(par.tsls, rdata, edata, par.pos)
    par <- par.tsls
    for(i in (ncut-1):0){
      v <- i/ncut * s0
      t2 <- try(par <- newton.raphson.LRT(par, rdata, edata, par.pos, v), silent = TRUE)
      if('try-error' %in% class(t2)){
        next
      }
    }

    if(!check.LRT.mle(par, rdata, edata, par.pos)){
      stop('Cannot find MLE')
    }
  }

  par

}
