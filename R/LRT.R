
LRT <- function(rdata, edata, par, se, c.lrt, level){

  tsls <- find.tsls(rdata, edata)
  par.pos <- align.parameter(tsls)$par.pos
  logL.obs.alt <- logL(par, rdata, edata, par.pos)
  rm(tsls, par.pos)

  logL.obs.null <- logL.null(rdata, edata)
  lrt.stat <- -2 * (logL.obs.null - logL.obs.alt)

  p.lrt <- pchisq(lrt.stat, df = 1, lower.tail = FALSE)
  ap.lrt <- pchisq(lrt.stat/c.lrt, df = 1, lower.tail = FALSE)

  ci <- LRT.conf.int(rdata, edata, par, se, c.lrt, level)

  list(p.lrt = p.lrt, ap.lrt = ap.lrt, c.lrt = c.lrt, ci = ci)

}

