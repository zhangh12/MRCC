
line.search <- function(par, rdata, edata, par.pos, delta){

  f0 <- neg.log.lik(par, rdata, edata, par.pos)
  f <- NULL
  for(k in 1:200){
    par1 <- par
    par1[] <- par - .01 * k * delta
    f <- c(f, neg.log.lik(par1, rdata, edata, par.pos))
  }

  min.f <- ifelse(all(is.na(f)), NA, min(f, na.rm = TRUE))
  if(!is.na(min.f) && min.f < f0){
    #message('line search updates')
    k <- which(!is.na(f) & f == min.f)[1]
    par[] <- par - .01 * k * delta
    status <- TRUE
  }else{
    #message('line search fails')
    status <- FALSE
  }

  list(par = par, status = status)

}