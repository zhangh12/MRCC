
neg.log.lik <- function(par, rdata, edata, par.pos){

  if(any(is.na(par))){
    return(NULL)
  }

  ne <- nrow(edata$data)
  nr <- nrow(rdata$data)

  a <- par['a']

  if('b' %in% par.pos$par.name){
    b <- par['b']
  }else{
    b <- NA
  }

  c <- par['c']
  alp.0 <- par['alp.0']

  if('alp.x' %in% par.pos$par.name){
    alp.x <- par[par.pos['alp.x', 'start']:par.pos['alp.x', 'end']]
  }else{
    alp.x <- NA
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

  bet.z <- par[par.pos['bet.z', 'start']:par.pos['bet.z', 'end']]

  lin <- a + alp.0 * bet.z + .5 * exp(c) * bet.z^2
  if(!is.na(b)){
    lin <- lin + b * bet.z
  }

  rg <- as.matrix(rdata$data[, rdata$vg, drop = FALSE])
  eg <- as.matrix(edata$data[, edata$vg, drop = FALSE])

  lin <- lin + rg %*% alp.g * bet.z
  if(length(rdata$vx) > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    lin <- lin + rx %*% (bet.x + bet.z * alp.x)
  }

  if(length(rdata$vy) > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
    lin <- lin + ry %*% bet.y
  }

  lin <- as.vector(lin)
  delta <- exp(lin)

  rd <- rdata$data[, rdata$vd]
  n1 <- sum(rd)
  n0 <- nr - n1

  res <- edata$data[, edata$vz] - alp.0 - eg %*% alp.g
  if(!is.na(b)){
    res <- res - (bet.z * exp(c) + b) * edata$data[, edata$vd]
  }

  if(length(edata$vx) > 0){
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    res <- res - ex %*% alp.x
  }

  -(sum(rd * lin) - nr * log(n0) - sum(log(1 + n1 / n0 * delta)) - ne * log(2*pi)/2 - ne * c/2 - exp(-c)/2 * sum(res^2))

}


