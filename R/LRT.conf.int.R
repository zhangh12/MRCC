
LRT.conf.int <- function(rdata, edata, par, se, c.lrt, level){

  tsls <- find.tsls(rdata, edata)
  par.pos <- align.parameter(tsls)$par.pos
  rm(tsls)

  logL.obs <- logL(par, rdata, edata, par.pos)
  crit <- logL.obs - 0.5 * c.lrt * qchisq(level, df = 1)

  name.bet.z <- paste0('bet.', edata$vz)
  se <- se[name.bet.z]
  b0 <- par[name.bet.z]
  names(b0) <- NULL
  names(se) <- NULL

  par0 <- par
  L0 <- logL(par0, rdata, edata, par.pos) - crit # L0 > 0
  if(L0 <= 0){
    stop('L0 <= 0')
  }

  x <- b0
  y <- L0
  L1 <- NULL
  for(i in c(0.5, 1:100)){
    b1 <- b0 - qnorm((1+level)/2) * se * i
    t1 <- try(par1 <- find.LRT.mle(rdata, edata, par.pos, b1))
    if('try-error' %in% class(t1)){
      next
    }
    L1 <- logL(par1, rdata, edata, par.pos) - crit
    x <- c(x, b1)
    y <- c(y, L1)
    if(L0 * L1 <= 0){
      break
    }
  }

  if(is.null(L1) | L0 * L1 > 0){
    lci <- -Inf
  }else{
    lci <- NULL
  }

  L2 <- NULL
  for(i in c(0.5, 1:100)){
    b2 <- b0 + qnorm((1+level)/2) * se * i
    t2 <- try(par2 <- find.LRT.mle(rdata, edata, par.pos, b2))
    if('try-error' %in% class(t2)){
      next
    }
    L2 <- logL(par2, rdata, edata, par.pos) - crit
    x <- c(x, b2)
    y <- c(y, L2)
    if(L0 * L2 <= 0){
      break
    }
  }

  if(is.null(L2) | L0 * L2 > 0){
    rci <- Inf
  }else{
    rci <- NULL
  }

  #plot(sort(x), y[order(x)],type='b',pch=20)
  #abline(h=0)

  # search for left end point of CI
  b00 <- b0
  if(is.null(lci)){
    while(TRUE){

      b <- (b0 + b1)/2
      parb <- find.LRT.mle(rdata, edata, par.pos, b)
      Lb <- logL(parb, rdata, edata, par.pos) - crit
      if(Lb > 1e-6){
        b0 <- b
        next
      }

      if(Lb < -1e-6){
        b1 <- b
        next
      }

      lci <- b
      break
    }
  }

  # search for right end point of CI
  b0 <- b00
  if(is.null(rci)){
    while(TRUE){

      b <- (b0 + b2)/2
      parb <- find.LRT.mle(rdata, edata, par.pos, b)
      Lb <- logL(parb, rdata, edata, par.pos) - crit
      if(Lb > 1e-6){
        b0 <- b
        next
      }

      if(Lb < -1e-6){
        b2 <- b
        next
      }

      rci <- b
      break
    }
  }

  c(LCL = lci, RCL = rci)

}