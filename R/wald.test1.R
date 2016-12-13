
wald.test1 <- function(rdata, edata, c.wald, level, mle.only = FALSE){

  tsls <- find.tsls(rdata, edata)

  ap <- align.parameter(tsls)
  par.tsls <- ap$par
  par.pos <- ap$par.pos

  t1 <- try(par <- find.mle.2S(par.tsls, rdata, edata, par.pos))
  if('try-error' %in% class(t1)){
    return(NULL)
  }

  max.logL <- logL(par, rdata, edata, par.pos)
  hess <- hessian(par, rdata, edata, par.pos)
  se1 <- sqrt(diag(solve(-hess)))

  wald.stat <- (par[paste0('bet.', edata$vz)]/se1[paste0('bet.', edata$vz)])^2

  if(mle.only){
    return(list(par = par, wald.stat = wald.stat, max.logL = max.logL))
  }

  hess0 <- hessian0(par, rdata, edata, par.pos)

  ev <- empirical.variance(par, rdata, edata, par.pos)
  tv <- theoretical.variance(par, rdata, edata, par.pos)
  se <- sqrt(diag(ev))
  se0 <- sqrt(diag(tv))

  se2 <- sqrt(diag(solve(-hess0)))

  summary.n <- data.frame(Estimate = par, SE = se1, SE0 = se0, stringsAsFactors = FALSE)

  rownames(summary.n) <- names(par)
  summary.n$z <- summary.n$Estimate / summary.n$SE
  summary.n$"Pr(>|z|)" <- pchisq(summary.n$z^2, df = 1, lower.tail = FALSE)

  ap.wald <- pchisq(wald.stat/c.wald, df = 1, lower.tail = FALSE)

  ci <- par[paste0('bet.', edata$vz)] + c(-1, 1) * se1[paste0('bet.', edata$vz)] * sqrt(c.wald * qchisq(level, df = 1))
  names(ci) <- c('LCL', 'RCL')

  gr <- score(par, rdata, edata, par.pos)

  #####

  list(par = par, se = se1, ap.wald = ap.wald, c.wald = c.wald, ci = ci, summary.n = summary.n, gr = gr, max.logL = max.logL)

}

