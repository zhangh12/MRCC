

mrcc <- function(rformula, rdata, eformula, edata, nb = 500, level = 0.95, start = NULL, b.ci = NULL, fig = FALSE){

  rdata <- parse.rdata(rformula, rdata)
  edata <- parse.edata(eformula, edata)

  exd <- extend.data(rdata, edata)

  rdata <- exd$rdata
  edata <- exd$edata
  summary.h <- exd$summary.h

  c.adj <- null.adjustment(rdata, edata, nb)
  res.tsls <- TSLS.test(rdata, edata, c.adj$c.tsls, level, FALSE, b.ci)
  res.wald <- wald.test(rdata, edata, c.adj$c.wald, level, FALSE, start, b.ci)
  if(is.null(res.wald)){
    return(NULL)
  }

  find.profile.mle(rdata, edata)

  #res.wald1 <- wald.test1(rdata, edata, c.adj$c.wald, level)
  res.lm <- LM.test(rdata, edata, res.tsls$par, res.tsls$se, level, b.ci, fig)
  res.lr <- LR.test(rdata, edata, res.wald$max.logL, res.lm$stat, res.tsls$par, res.tsls$se, c.adj$c.lrt, level, b.ci, fig)

  res.wald <- rescale.wald.ci(res.wald, res.lm, edata, level, b.ci)

  #return(list(res.wald = res.wald, res.tsls = res.tsls))
  return(list(res.wald = res.wald, res.lr = res.lr, res.lm = res.lm, res.tsls = res.tsls))

}

