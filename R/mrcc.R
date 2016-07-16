

mrcc <- function(rformula, rdata, eformula, edata){

  rdata <- parse.rdata(rformula, rdata)
  edata <- parse.edata(eformula, edata)

  exd <- extend.data(rdata, edata)

  rdata <- exd$rdata
  edata <- exd$edata

  mrcc.test(rdata, edata)

}

