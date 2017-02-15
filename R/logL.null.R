
# calculate max log-likelihood under the null (bet = 0)
logL.null <- function(rdata, edata){

  null.tsls <- find.null.tsls(rdata, edata)
  ap <- align.parameter(null.tsls)
  par.null <- ap$par
  par.pos <- ap$par.pos

  logL(par.null, rdata, edata, par.pos)

}

