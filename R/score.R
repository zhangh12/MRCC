

score <- function(par, rdata, edata, par.pos){

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
    nx <- length(alp.x)
  }else{
    alp.x <- NA
    nx <- 0
  }

  alp.g <- par[par.pos['alp.g', 'start']:par.pos['alp.g', 'end']]
  ng <- length(alp.g)

  if('bet.y' %in% par.pos$par.name){
    bet.y <- par[par.pos['bet.y', 'start']:par.pos['bet.y', 'end']]
    ny <- length(bet.y)
  }else{
    bet.y <- NA
    ny <- 0
  }

  if('bet.x' %in% par.pos$par.name){
    bet.x <- par[par.pos['bet.x', 'start']:par.pos['bet.x', 'end']]
  }else{
    bet.x <- NA
  }

  bet.z <- par[par.pos['bet.z', 'start']:par.pos['bet.z', 'end']]

  ##################

  rg <- as.matrix(rdata$data[, rdata$vg, drop = FALSE])
  eg <- as.matrix(edata$data[, edata$vg, drop = FALSE])
  eta.g <- as.vector(rg %*% alp.g)
  eta <- eta.g

  if(nx > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    eta.x <- as.vector(rx %*% alp.x)
    eta <- eta + eta.x
  }
  eta <- as.vector(eta)

  if(ny > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
  }

  rd <- as.vector(rdata$data[, rdata$vd])
  n1 <- sum(rd)
  n0 <- nr - n1

  lin <- a

  lin <- lin + rg %*% alp.g * bet.z
  if(nx > 0){
    lin <- lin + rx %*% (bet.x + bet.z * alp.x)
  }

  if(ny > 0){
    lin <- lin + ry %*% bet.y
  }

  lin <- as.vector(lin)

  delta <- exp(lin)

  p <- 1/ n0 / (1 + n1/n0 * delta) # 1 / (n0 + n1 * delta)
  Delta <- -1 + 1 / (1 + n1/n0 * delta) # -n1 p delta
  xi <- Delta * (1 + Delta)

  res <- edata$data[, edata$vz] - alp.0 - eg %*% alp.g
  if(!is.na(b)){
    res <- res - (bet.z * exp(c) + b) * edata$data[, edata$vd]
  }

  if(nx > 0){
    res <- res - ex %*% alp.x
  }
  res <- as.vector(res)

  one <- rep(1, ne)

  name.alp.g <- paste0('alp.', rdata$vg)
  name.alp.x <- paste0('alp.', rdata$vx)
  name.bet.x <- paste0('bet.', rdata$vx)
  name.bet.y <- paste0('bet.', rdata$vy)
  name.bet.z <- paste0('bet.', edata$vz)

  ############################
  ## calculate score vector ##
  ############################

  gr <- rep(NA, length(par))
  names(gr) <- names(par)

  gr['c'] <- -ne/2 + exp(-c)/2 * sum(res^2)

  gr['alp.0'] <- exp(-c) * sum(res)

  if(nx > 0){
    gr[name.alp.x] <- as.vector(bet.z * (t(rx) %*% (rd + Delta)) + exp(-c) * (t(ex) %*% res))
  }

  gr[name.alp.g] <- as.vector(bet.z * (t(rg) %*% (rd + Delta)) + exp(-c) * (t(eg) %*% res))

  gr['a'] <- sum(rd + Delta)

  if(nx > 0){
    gr[name.bet.x] <- as.vector(t(rx) %*% (rd + Delta))
  }

  if(ny > 0){
    gr[name.bet.y] <- as.vector(t(ry) %*% (rd + Delta))
  }

  gr[name.bet.z] <- sum((rd + Delta) * eta)

  gr <- - gr

  gr


}