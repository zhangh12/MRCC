
logL.LRT <- function(par, rdata, edata, par.pos, bet){

  if(any(is.na(par))){
    return(NULL)
  }

  ne <- nrow(edata$data)
  nr <- nrow(rdata$data)

  a <- par[par.pos['a', 'start']]

  if('b' %in% names(par)){
    b <- par[par.pos['b', 'start']]
    est.b <- TRUE
  }else{
    b <- 0
    est.b <- FALSE
  }

  c <- par[par.pos['c', 'start']]
  alp.0 <- par[par.pos['alp.0', 'start']]

  if('alp.x' %in% par.pos$par.name){
    alp.x <- par[par.pos['alp.x', 'start']:par.pos['alp.x', 'end']]
  }else{
    alp.x <- NA
  }

  if('alp.h' %in% par.pos$par.name){
    alp.h <- par[par.pos['alp.h', 'start']:par.pos['alp.h', 'end']]
  }else{
    alp.h <- NA
  }

  alp.g <- par[par.pos['alp.g', 'start']:par.pos['alp.g', 'end']]

  if('bet.y' %in% par.pos$par.name){
    bet.y <- par[par.pos['bet.y', 'start']:par.pos['bet.y', 'end']]
  }else{
    bet.y <- NA
  }

  if('bet.x' %in% par.pos$par.name){
    bet.x <- par[par.pos['bet.x', 'start']:par.pos['bet.x', 'end']]
  }else{
    bet.x <- NA
  }

  bet.z <- bet

  lin <- a

  rg <- as.matrix(rdata$data[, rdata$vg, drop = FALSE])
  eg <- as.matrix(edata$data[, edata$vg, drop = FALSE])

  lin <- lin + rg %*% alp.g * bet.z
  if(length(rdata$vx) > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    lin <- lin + rx %*% bet.x
  }

  if(length(rdata$vy) > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
    lin <- lin + ry %*% bet.y
  }

  lin <- as.vector(lin)
  delta <- exp(lin)

  rd <- as.vector(rdata$data[, rdata$vd])
  n1 <- sum(rd)
  n0 <- nr - n1

  r <- edata$data[, edata$vz] - alp.0 - eg %*% alp.g - (bet.z * c + b) * edata$data[, edata$vd]

  if(length(edata$vx) > 0){
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    r <- r - ex %*% alp.x
  }

  if(length(edata$vh) > 0){
    eh <- as.matrix(edata$data[, edata$vh, drop = FALSE])
    r <- r - eh %*% alp.h
  }
  r <- as.vector(r)

  l <- sum(rd * lin) - sum(log(n0 + n1 * delta)) - ne * log(2*pi)/2 - ne * log(c)/2 - 1/c/2 * sum(r^2)
  names(l) <- NULL

  l

}


