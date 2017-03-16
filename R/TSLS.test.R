
# wald's test for two-stage method
TSLS.test <- function(rdata, edata, omega, level, b.ci = NULL){

  tsls <- find.tsls(rdata, edata, TRUE)

  ap <- align.parameter(tsls)
  par <- ap$par
  par.pos <- ap$par.pos

  sc <- score.TSLS(par, rdata, edata, par.pos) # close to zero
  #hess <- hessian.TSLS(par, rdata, edata, par.pos)
  hess0 <- hessian0.TSLS(par, rdata, edata, omega, par.pos)
  #info <- fisher.info.TSLS(par, rdata, edata, omega, par.pos)
  info0 <- fisher.info0.TSLS(par, rdata, edata, omega, par.pos)

  cov <- solve(hess0) %*% info0 %*% t(solve(hess0))
  se <- sqrt(diag(cov))
  names(se) <- names(par)

  name.bet.z <- paste0('bet.', edata$vz)
  tsls.stat <- (par[name.bet.z]/se[name.bet.z])^2

  if('b' %in% names(par)){
    par['b'] <- par['b'] - par[name.bet.z] * par['c']
    v <- c(-par[name.bet.z], 1, -par['c'])
    se['b'] <- sqrt(t(v) %*% cov[c('c', 'b', name.bet.z), c('c', 'b', name.bet.z)] %*% v)
  }

  summary <- data.frame(Estimate = par, SE = se, stringsAsFactors = FALSE)

  rownames(summary) <- names(par)
  summary$z <- summary$Estimate / summary$SE
  summary$"Pr(>|z|)" <- pchisq(summary$z^2, df = 1, lower.tail = FALSE)

  p.tsls <- pchisq(tsls.stat, df = 1, lower.tail = FALSE)

  ci <- par[name.bet.z] + c(-1, 1) * se[name.bet.z] * sqrt(qchisq(level, df = 1))
  names(ci) <- c('LCL', 'RCL')

  in.ci <- ifelse(is.null(b.ci), NA, ci['LCL'] <= b.ci && ci['RCL'] >= b.ci)

  list(par = par, se = se, p.tsls = p.tsls, gr = sc,
       ci = ci, summary = summary, in.ci = in.ci)

}
