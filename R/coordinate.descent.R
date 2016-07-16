
coordinate.descent <- function(par, rdata, edata, par.pos){

  npar <- length(par)

  f0 <- neg.log.lik(par, rdata, edata, par.pos)
  sc <- score(par, rdata, edata, par.pos)
  f <- NULL
  optim.step <- NULL
  for(i in 1:npar){
    g <- NULL
    s <- seq(0, 1, by = 0.02)
    s <- s[s != 0]
    for(j in s){
      par1 <- par
      par1[i] <- par[i] - j * sc[i]
      g <- c(g, neg.log.lik(par1, rdata, edata, par.pos))
    }

    g[is.infinite(g)] <- Inf
    min.g <- ifelse(all(is.na(g)), Inf, min(g, na.rm = TRUE))
    optim.step <- c(optim.step, ifelse(is.infinite(min.g), NA, s[which(!is.na(g) & g == min.g)]))

    f <- c(f, min.g)
  }

  min.f <- ifelse(all(is.na(f)), NA, min(f, na.rm = TRUE))
  if(!is.na(min.f) && min.f < f0){
    #message('coordinate-wise updates')
    k <- which(!is.na(f) & f == min.f)
    par[k] <- par[k] - optim.step[k] * sc[k]
    status <- TRUE
  }else{
    #message('coordinate descent fails')
    status <- FALSE
  }

  list(par = par, status = status)

}