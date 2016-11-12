

score <- function(par, rdata, edata, par.pos, v = 0){

  if(any(is.na(par))){
    return(NULL)
  }

  ne <- nrow(edata$data)
  nr <- nrow(rdata$data)

  a <- par['a']

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
  eta <- as.vector(rg %*% alp.g)

  if(nx > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    #eta.x <- as.vector(rx %*% alp.x)
  }

  if(ny > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
  }

  rd <- as.vector(rdata$data[, rdata$vd])
  n1 <- sum(rd)
  n0 <- nr - n1

  lin <- a + rg %*% alp.g * bet.z
  if(nx > 0){
    lin <- lin + rx %*% bet.x
  }

  if(ny > 0){
    lin <- lin + ry %*% bet.y
  }

  lin <- as.vector(lin)

  delta <- exp(lin)

  p <- 1/ n0 / (1 + n1/n0 * delta) # 1 / (n0 + n1 * delta)
  Delta <- -1 + 1 / (1 + n1/n0 * delta) # -n1 p delta
  xi <- Delta * (1 + Delta)

  r <- edata$data[, edata$vz] - alp.0 - eg %*% alp.g

  if(nx > 0){
    r <- r - ex %*% alp.x
  }
  r <- as.vector(r)

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

  gr['c'] <- -ne/c/2 + 1/c^2/2 * sum(r^2)

  gr['alp.0'] <- 1/c * sum(r)

  if(nx > 0){
    gr[name.alp.x] <- as.vector(1/c * t(ex) %*% r)
  }

  gr[name.alp.g] <- as.vector(bet.z * (t(rg) %*% (rd + Delta)) + 1/c * (t(eg) %*% r))

  gr['a'] <- sum(rd + Delta)

  if(nx > 0){
    gr[name.bet.x] <- as.vector(t(rx) %*% (rd + Delta))
  }

  if(ny > 0){
    gr[name.bet.y] <- as.vector(t(ry) %*% (rd + Delta))
  }

  gr[name.bet.z] <- sum((rd + Delta) * eta)

  gr - v


}