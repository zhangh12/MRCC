
LR.test <- function(rdata, edata, logL.obs.alt, lm.stat, par, se, level, b.ci, fig){

  logL.obs.null <- logL.null(rdata, edata)
  lr.stat <- -2 * (logL.obs.null - logL.obs.alt)

  c.lrt <- lr.stat/lm.stat

  p.lrt <- pchisq(lr.stat/c.lrt, df = 1, lower.tail = FALSE)

  ci <- LRT.conf.int(rdata, edata, logL.obs.alt, par, se, c.lrt, level, fig)

  in.ci <- check.LRT.ci(rdata, edata, logL.obs.alt, c.lrt, level, b.ci)

  list(p.lrt = p.lrt, c.lrt = c.lrt, ci = ci, in.ci = in.ci, stat = lr.stat)

}

