

mrcc <- function(rformula, rdata, eformula, edata, nb = 500, level = 0.95){

  rdata <- parse.rdata(rformula, rdata)
  edata <- parse.edata(eformula, edata)

  exd <- extend.data(rdata, edata)

  rdata <- exd$rdata
  edata <- exd$edata
  summary.h <- exd$summary.h

  c.adj <- null.adjustment(rdata, edata, nb)
  res.wald <- wald.test(rdata, edata, c.adj$c.wald, level)
  if(is.null(res.wald)){
    return(NULL)
  }
  #res.wald1 <- wald.test1(rdata, edata, c.adj$c.wald, level)
  res.lrt <- LRT(rdata, edata, res.wald$par, res.wald$se, c.adj$c.lrt, level)
  res.tsls <- TSLS.test(rdata, edata, c.adj$c.tsls, level)

  summary.mrcc(rdata, edata, c.adj, res.wald, res.lrt, res.tsls)

  list(summary.h = summary.h, res.wald = res.wald, res.lrt = res.lrt, res.tsls = res.tsls, c.adj = c.adj)

}

