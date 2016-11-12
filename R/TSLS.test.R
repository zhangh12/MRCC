
TSLS.test <- function(rdata, edata){

  tsls <- find.tsls(rdata, edata)

  ap <- align.parameter(tsls)
  par <- ap$par
  par.pos <- ap$par.pos

  hess <- hessian0.TSLS(par, rdata, edata, par.pos)
  info <- fisher.info0.TSLS(par, rdata, edata, par.pos)

  cov <- solve(hess) %*% info %*% t(solve(hess))
  se <- sqrt(diag(cov))
  names(se) <- names(par)

  summary.n <- data.frame(Estimate = par, SE = se, stringsAsFactors = FALSE)

  rownames(summary.n) <- names(par)
  summary.n$z <- summary.n$Estimate / summary.n$SE
  summary.n$"Pr(>|z|)" <- pchisq(summary.n$z^2, df = 1, lower.tail = FALSE)

  summary.n

}
