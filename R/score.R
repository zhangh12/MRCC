

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

  ##################

  rg <- as.matrix(rdata$data[, rdata$vg, drop = FALSE])
  eg <- as.matrix(edata$data[, edata$vg, drop = FALSE])
  prd <- as.vector(alp.0 + exp(c) * bet.z + rg %*% alp.g)

  if(length(rdata$vx) > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    prd <- prd + as.vector(rx %*% alp.x)
  }
  prd <- as.vector(prd)

  if(length(rdata$vy) > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
  }

  rd <- as.vector(rdata$data[, rdata$vd])
  n1 <- sum(rd)
  n0 <- nr - n1

  lin <- a + alp.0 * bet.z + .5 * exp(c) * bet.z^2
  if(!is.na(b)){
    lin <- lin + b * bet.z
  }

  lin <- lin + rg %*% alp.g * bet.z
  if(length(rdata$vx) > 0){
    lin <- lin + rx %*% (bet.x + bet.z * alp.x)
  }

  if(length(rdata$vy) > 0){
    lin <- lin + ry %*% bet.y
  }

  lin <- as.vector(lin)

  delta <- exp(lin)

  p <- 1/ n0 / (1 + n1/n0 * delta) # 1 / (n0 + n1 * delta)
  d.hat <- 1 - 1 / (1 + n1/n0 * delta) # n1 * delta / (n0 + n1 * delta) = 1 - n0 / (n0 + n1 * delta)
  Delta <- rd - d.hat
  xi <- -d.hat * (1 - d.hat)

  res <- edata$data[, edata$vz] - alp.0 - eg %*% alp.g
  if(!is.na(b)){
    res <- res - (bet.z * exp(c) + b) * edata$data[, edata$vd]
  }

  if(length(edata$vx) > 0){
    res <- res - ex %*% alp.x
  }
  res <- as.vector(res)

  ############################
  ## calculate score vector ##
  ############################

  gr <- rep(NA, length(par))
  names(gr) <- names(par)

  ell.a <- sum(Delta)
  gr['a'] <- ell.a

  ell.c <- .5 * exp(c) * bet.z^2 * sum(Delta) - ne/2 + exp(-c)/2 * sum(res^2)
  gr['c'] <- ell.c

  if(length(rdata$vx) > 0){
    ell.bet.x <- t(rx) %*% Delta
    ell.bet.x <- as.vector(ell.bet.x)
    gr[paste0('bet.', rdata$vx)] <- ell.bet.x
  }

  if(length(rdata$vy) > 0){
    ell.bet.y <- t(ry) %*% Delta
    ell.bet.y <- as.vector(ell.bet.y)
    gr[paste0('bet.', rdata$vy)] <- ell.bet.y
  }

  ell.bet.z <- sum(Delta * prd)
  gr[paste0('bet.', edata$vz)] <- ell.bet.z

  ell.alp.0 <- bet.z * sum(Delta) + exp(-c) * sum(res)
  gr['alp.0'] <- ell.alp.0

  if(length(rdata$vx) > 0){
    ell.alp.x <- bet.z * t(rx) %*% Delta + exp(-c) * t(ex) %*% res
    ell.alp.x <- as.vector(ell.alp.x)
    gr[paste0('alp.', edata$vx)] <- ell.alp.x
  }

  ell.alp.g <- bet.z * t(rg) %*% Delta + exp(-c) * t(eg) %*% res
  ell.alp.g <- as.vector(ell.alp.g)
  gr[paste0('alp.', edata$vg)] <- ell.alp.g

  gr <- - gr

  gr


}