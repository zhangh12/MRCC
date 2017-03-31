
wald.test <- function(rdata, edata, omega, level, start = NULL, b.ci = NULL){

  if(is.null(start)){
    tsls <- find.tsls(rdata, edata)
  }else{
    tsls <- start
  }

  ap <- align.parameter(tsls)
  par.tsls <- ap$par
  par.pos <- ap$par.pos

  t0 <- try(par <- find.mle(par.tsls, rdata, edata, par.pos))
  if('try-error' %in% class(t0)){
    return(NULL)
  }

  #logL.ind(par, rdata, edata, par.pos)
  #tau.a33(par, rdata, edata, par.pos)

  name.bet.z <- paste0('bet.', edata$vz)
  max.logL <- logL(par, rdata, edata, par.pos)

  se <- working.variance(par, rdata, edata, omega, par.pos)

  emp.se <- empirical.variance(par, rdata, edata, omega, par.pos)
  se0 <- theoretical.variance(par, rdata, edata, omega, par.pos)

  summary <- data.frame(Estimate = par, SE = se, stringsAsFactors = FALSE)

  rownames(summary) <- names(par)
  summary$z <- summary$Estimate / summary$SE
  summary$"Pr(>|z|)" <- pchisq(summary$z^2, df = 1, lower.tail = FALSE)
  wald.stat <- (par[name.bet.z]/se[name.bet.z])^2
  p.wald <- pchisq(wald.stat, df = 1, lower.tail = FALSE)

  ci <- par[name.bet.z] + c(-1, 1) * se[name.bet.z] * sqrt(qchisq(level, df = 1))
  names(ci) <- c('LCL', 'RCL')

  gr <- score(par, rdata, edata, par.pos)

  #####

  in.ci <- ifelse(is.null(b.ci), NA, ci['LCL'] <= b.ci && ci['RCL'] >= b.ci)

  list(par = par, se = se, p.wald = p.wald, ci = ci,
       summary = summary, in.ci = in.ci, gr = gr, max.logL = max.logL,
       stat = wald.stat)

}

