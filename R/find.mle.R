
find.mle <- function(par.tsls, rdata, edata, par.pos, ncut = 50){

  t1 <- try(par <- newton.raphson(par.tsls, rdata, edata, par.pos), silent = TRUE)
  if(('try-error' %in% class(t1)) | !check.mle(par, rdata, edata, par.pos)){
    s0 <- score(par.tsls, rdata, edata, par.pos)
    par <- par.tsls
    for(i in (ncut-1):0){
      v <- i/ncut * s0
      t2 <- try(par <- newton.raphson(par, rdata, edata, par.pos, v), silent = TRUE)
      if('try-error' %in% class(t2)){
        next
      }
    }

    if(!check.mle(par, rdata, edata, par.pos)){
      stop('Cannot find MLE')
    }
  }

  par

}