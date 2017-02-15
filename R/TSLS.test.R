
# perform wald's test for two-stage method
# bootstrap adjusted if c.tsls is provided
TSLS.test <- function(rdata, edata, c.tsls, level, tsls.only = FALSE, b.ci = NULL){

  tsls <- find.tsls(rdata, edata)

  ap <- align.parameter(tsls)
  par <- ap$par
  par.pos <- ap$par.pos

  hess <- hessian0.TSLS(par, rdata, edata, par.pos)
  info <- fisher.info0.TSLS(par, rdata, edata, par.pos)

  cov <- solve(hess) %*% info %*% t(solve(hess))
  se <- sqrt(diag(cov))
  names(se) <- names(par)

  tsls.stat <- (par[paste0('bet.', edata$vz)]/se[paste0('bet.', edata$vz)])^2

  if(tsls.only){
    return(list(tsls.stat = tsls.stat))
  }

  summary.n <- data.frame(Estimate = par, SE = se, stringsAsFactors = FALSE)

  rownames(summary.n) <- names(par)
  summary.n$z <- summary.n$Estimate / summary.n$SE
  summary.n$"Pr(>|z|)" <- pchisq(summary.n$z^2, df = 1, lower.tail = FALSE)

  p.tsls <- pchisq(tsls.stat, df = 1, lower.tail = FALSE)
  ap.tsls <- pchisq(tsls.stat/c.tsls, df = 1, lower.tail = FALSE)

  ci <- par[paste0('bet.', edata$vz)] + c(-1, 1) * se[paste0('bet.', edata$vz)] * sqrt(c.tsls * qchisq(level, df = 1))
  names(ci) <- c('LCL', 'RCL')

  in.ci <- ifelse(is.null(b.ci), NA, ci['LCL'] <= b.ci && ci['RCL'] >= b.ci)

  list(par = par, se = se, p.tsls = p.tsls, ap.tsls = ap.tsls,
       c.tsls = c.tsls, ci = ci, summary.n = summary.n, in.ci = in.ci)

}
