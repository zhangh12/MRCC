

hessian.LRT <- function(par, rdata, edata, par.pos){

  h <- hessian(par, rdata, edata, par.pos)

  name.bet.z <- paste0('bet.', edata$vz)

  id <- which(colnames(h) == name.bet.z)
  h <- h[-id, -id]

  h

}

