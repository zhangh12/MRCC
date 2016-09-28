

mrcc <- function(rformula, rdata, eformula, edata, tpar = NULL){

  rdata <- parse.rdata(rformula, rdata)
  edata <- parse.edata(eformula, edata)

  exd <- extend.data(rdata, edata)

  rdata <- exd$rdata
  edata <- exd$edata

  res <- mrcc.test(rdata, edata, tpar)
  res$TSLS <- TSLS.test(rdata, edata, tpar)

  res

}

