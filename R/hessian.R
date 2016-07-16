

hessian <- function(par, rdata, edata, par.pos){

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

  one <- rep(1, ne)

  ###########################
  ## calculate Hess matrix ##
  ###########################

  hess <- matrix(NA, nrow = length(par), ncol = length(par))
  rownames(hess) <- names(par)
  colnames(hess) <- names(par)

  ell.a.a <- sum(xi)
  hess['a', 'a'] <- ell.a.a

  ell.a.c <- exp(c)/2 * bet.z^2 * sum(xi)
  hess['a', 'c'] <- ell.a.c

  if(length(rdata$vx) > 0){
    ell.a.bet.x <- t(t(rx) %*% xi)
    hess['a', paste0('bet.', rdata$vx)] <- ell.a.bet.x
  }

  if(length(rdata$vy) > 0){
    ell.a.bet.y <- t(t(ry) %*% xi)
    hess['a', paste0('bet.', rdata$vy)] <- ell.a.bet.y
  }

  ell.a.bet.z <- sum(xi * prd)
  hess['a', paste0('bet.', edata$vz)] <- ell.a.bet.z

  ell.a.alp.0 <- bet.z * sum(xi)
  hess['a', 'alp.0'] <- ell.a.alp.0

  if(length(rdata$vx) > 0){
    ell.a.alp.x <- bet.z * t(t(rx) %*% xi)
    hess['a', paste0('alp.', rdata$vx)] <- ell.a.alp.x
  }

  ell.a.alp.g <- bet.z * t(t(rg) %*% xi)
  hess['a', paste0('alp.', rdata$vg)] <- ell.a.alp.g

  #########

  ell.c.c <- exp(c)/2 * bet.z^2 * sum(Delta) + exp(2*c)/4 * bet.z^4 * sum(xi) - exp(-c)/2 * sum(res^2)
  hess['c', 'c'] <- ell.c.c

  if(length(rdata$vx) > 0){
    ell.c.bet.x <- exp(c)/2 * bet.z^2 * t(t(rx) %*% xi)
    hess['c', paste0('bet.', rdata$vx)] <- ell.c.bet.x
  }

  if(length(rdata$vy) > 0){
    ell.c.bet.y <- exp(c)/2 * bet.z^2 * t(t(ry) %*% xi)
    hess['c', paste0('bet.', rdata$vy)] <- ell.c.bet.y
  }

  ell.c.bet.z <- exp(c) * bet.z * sum(Delta) + exp(c)/2 * bet.z^2 * sum(xi * prd)
  hess['c', paste0('bet.', edata$vz)] <- ell.c.bet.z

  ell.c.alp.0 <- exp(c)/2 * bet.z^3 * sum(xi) - exp(-c) * sum(res)
  hess['c', 'alp.0'] <- ell.c.alp.0

  if(length(rdata$vx) > 0){
    ell.c.alp.x <- exp(c)/2 * bet.z^3 * t(t(rx) %*% xi) - exp(-c) * t(t(ex) %*% res)
    hess['c', paste0('alp.', rdata$vx)] <- ell.c.alp.x
  }

  ell.c.alp.g <- exp(c)/2 * bet.z^3 * t(t(rg) %*% xi) - exp(-c) * t(t(eg) %*% res)
  hess['c', paste0('alp.', rdata$vg)] <- ell.c.alp.g

  ############

  if(length(rdata$vx) > 0){
    ell.bet.x.bet.x <- t(rx) %*% (xi * rx)
    hess[paste0('bet.', rdata$vx), paste0('bet.', rdata$vx)] <- ell.bet.x.bet.x

    if(length(rdata$vy) > 0){
      ell.bet.x.bet.y <- t(rx) %*% (xi * ry)
      hess[paste0('bet.', rdata$vx), paste0('bet.', rdata$vy)] <- ell.bet.x.bet.y
    }

    ell.bet.x.bet.z <- t(rx) %*% (xi * prd)
    hess[paste0('bet.', rdata$vx), paste0('bet.', edata$vz)] <- ell.bet.x.bet.z

    ell.bet.x.alp.0 <- bet.z * t(rx) %*% xi
    hess[paste0('bet.', rdata$vx), 'alp.0'] <- ell.bet.x.alp.0

    ell.bet.x.alp.x <- bet.z * t(rx) %*% (xi * rx)
    hess[paste0('bet.', rdata$vx), paste0('alp.', edata$vx)] <- ell.bet.x.alp.x

    ell.bet.x.alp.g <- bet.z * t(rx) %*% (xi * rg)
    hess[paste0('bet.', rdata$vx), paste0('alp.', rdata$vg)] <- ell.bet.x.alp.g
  }

  #############

  if(length(rdata$vy) > 0){
    ell.bet.y.bet.y <- t(ry) %*% (xi * ry)
    hess[paste0('bet.', rdata$vy), paste0('bet.', rdata$vy)] <- ell.bet.y.bet.y

    ell.bet.y.bet.z <- t(ry) %*% (xi * prd)
    hess[paste0('bet.', rdata$vy), paste0('bet.', edata$vz)] <- ell.bet.y.bet.z

    ell.bet.y.alp.0 <- bet.z * t(ry) %*% xi
    hess[paste0('bet.', rdata$vy), 'alp.0'] <- ell.bet.y.alp.0

    if(length(rdata$vx) > 0){
      ell.bet.y.alp.x <- bet.z * t(ry) %*% (xi * rx)
      hess[paste0('bet.', rdata$vy), paste0('alp.', rdata$vx)] <- ell.bet.y.alp.x
    }

    ell.bet.y.alp.g <- bet.z * t(ry) %*% (xi * rg)
    hess[paste0('bet.', rdata$vy), paste0('alp.', rdata$vg)] <- ell.bet.y.alp.g
  }

  #############

  ell.bet.z.bet.z <- exp(c) * sum(Delta) + sum(xi * prd^2)
  hess[paste0('bet.', edata$vz), paste0('bet.', edata$vz)] <- ell.bet.z.bet.z

  ell.bet.z.alp.0 <- sum(Delta) + bet.z * sum(xi * prd)
  hess[paste0('bet.', edata$vz), 'alp.0'] <- ell.bet.z.alp.0

  if(length(rdata$vx) > 0){
    ell.bet.z.alp.x <- t(t(rx) %*% Delta) + bet.z * t(t(rx) %*% (xi * prd))
    hess[paste0('bet.', edata$vz), paste0('alp.', rdata$vx)] <- ell.bet.z.alp.x
  }

  ell.bet.z.alp.g <- t(t(rg) %*% Delta) + bet.z * t(t(rg) %*% (xi * prd))
  hess[paste0('bet.', edata$vz), paste0('alp.', rdata$vg)] <- ell.bet.z.alp.g

  #############

  ell.alp.0.alp.0 <- bet.z^2 * sum(xi) - ne * exp(-c)
  hess['alp.0', 'alp.0'] <- ell.alp.0.alp.0

  if(length(rdata$vx) > 0){
    ell.alp.0.alp.x <- bet.z^2 * t(t(rx) %*% xi) - exp(-c) * t(t(ex) %*% one)
    hess['alp.0', paste0('alp.', rdata$vx)] <- ell.alp.0.alp.x
  }

  ell.alp.0.alp.g <- bet.z^2 * t(t(rg) %*% xi) - exp(-c) * t(t(eg) %*% one)
  hess['alp.0', paste0('alp.', rdata$vg)] <- ell.alp.0.alp.g

  #############

  if(length(rdata$vx) > 0){
    ell.alp.x.alp.x <- bet.z^2 * t(rx) %*% (xi * rx) - exp(-c) * t(ex) %*% ex
    hess[paste0('alp.', rdata$vx), paste0('alp.', rdata$vx)] <- ell.alp.x.alp.x

    ell.alp.x.alp.g <- bet.z^2 * t(rx) %*% (xi * rg) - exp(-c) * t(ex) %*% eg
    hess[paste0('alp.', rdata$vx), paste0('alp.', rdata$vg)] <- ell.alp.x.alp.g
  }

  #############

  ell.alp.g.alp.g <- bet.z^2 * t(rg) %*% (xi * rg) - exp(-c) * t(eg) %*% eg
  hess[paste0('alp.', rdata$vg), paste0('alp.', rdata$vg)] <- ell.alp.g.alp.g

  #############

  for(i in 1:nrow(hess)){
    for(j in 1:ncol(hess)){
      if(is.na(hess[i, j])){
        hess[i, j] <- hess[j, i]
      }
    }
  }

  hess <- -hess

  hess <- (hess + t(hess))/2

  hess

}

