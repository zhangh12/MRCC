

mrcc <- function(rformula, rdata, eformula, edata, level = 0.95, ce = FALSE, start = NULL, b.ci = NULL, fig = FALSE){

  rdata <- parse.rdata(rformula, rdata)
  edata <- parse.edata(eformula, edata)

  exd <- extend.data(rdata, edata)

  rdata <- exd$rdata
  edata <- exd$edata
  omega <- exd$omega

  res.tsls <- TSLS.test(rdata, edata, omega, level, b.ci)
  res.wald <- wald.test(rdata, edata, omega, level, start, b.ci)
  if(is.null(res.wald)){
    return(NULL)
  }

  #find.profile.mle(rdata, edata)

  res.lm <- LM.test(rdata, edata, omega, res.tsls$par, res.tsls$se, level, b.ci, fig)
  res.lr <- LR.test(rdata, edata, res.wald$max.logL, res.lm$stat, res.tsls$par, res.tsls$se, level, b.ci, fig)

  #res.wald <- rescale.wald.ci(res.wald, res.lm, edata, level, b.ci)

  list(wald = res.wald, lr = res.lr, lm = res.lm, tsls = res.tsls, version = packageVersion('MRCC'))

}

