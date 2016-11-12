

score.LRT <- function(par, rdata, edata, par.pos, v = 0){

  sc <- score(par, rdata, edata, par.pos)

  name.bet.z <- paste0('bet.', edata$vz)

  id <- which(names(sc) == name.bet.z)
  sc <- sc[-id] - v
  sc

}
