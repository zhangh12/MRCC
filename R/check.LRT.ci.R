
check.LRT.ci <- function(rdata, edata, logL.obs.alt, c.lrt, level, b.ci){

  if(is.null(b.ci)){
    return(NA)
  }

  par.pos <- align.parameter(find.tsls(rdata, edata))$par.pos
  crit <- c.lrt * qchisq(level, df = 1)
  name.bet.z <- paste0('bet.', edata$vz)

  if(abs(b.ci) < 1e-3){
    par.mle <- align.parameter(find.null.tsls(rdata, edata))$par
  }else{
    par.mle <-  try(find.LRT.mle(rdata, edata, b.ci), silent = TRUE)
    if('try-error' %in% class(par.mle)){
      return(NA)
    }
    par.mle <- c(par.mle, b.ci)
    names(par.mle)[par.pos[name.bet.z, 'start']] <- name.bet.z
  }

  lr.stat <- 2 * (logL.obs.alt - logL(par.mle, rdata, edata, par.pos))

  ifelse(lr.stat <= crit, TRUE, FALSE)

}

