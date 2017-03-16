

hessian.TSLS <- function(par, rdata, edata, par.pos){

  if(any(is.na(par))){
    return(NULL)
  }

  ne <- nrow(edata$data)
  nr <- nrow(rdata$data)

  a <- par['a']

  if('b' %in% names(par)){
    b <- par['b']
    est.b <- TRUE
  }else{
    b <- 0
    est.b <- FALSE
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

  if('alp.h' %in% par.pos$par.name){
    alp.h <- par[par.pos['alp.h', 'start']:par.pos['alp.h', 'end']]
    nh <- length(alp.h)
  }else{
    alp.h <- NA
    nh <- 0
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

  if(nh > 0){
    eh <- as.matrix(edata$data[, edata$vh, drop = FALSE])
  }

  if(ny > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
  }

  rd <- as.vector(rdata$data[, rdata$vd])
  n1 <- sum(rd)
  n0 <- nr - n1

  ed <- as.vector(edata$data[, edata$vd])
  ez <- as.vector(edata$data[, edata$vz])

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

  r <- ez - alp.0 - eg %*% alp.g - b * ed

  if(nx > 0){
    r <- r - ex %*% alp.x
  }

  if(nh > 0){
    r <- r - eh %*% alp.h
  }
  r <- as.vector(r)

  one <- rep(1, ne)

  name.alp.g <- paste0('alp.', edata$vg)
  name.alp.x <- paste0('alp.', edata$vx)
  name.alp.h <- paste0('alp.', edata$vh)
  name.bet.x <- paste0('bet.', rdata$vx)
  name.bet.y <- paste0('bet.', rdata$vy)
  name.bet.z <- paste0('bet.', edata$vz)


  ###########################
  ## calculate Hess matrix ##
  ###########################

  hess <- matrix(0, nrow = length(par), ncol = length(par))
  rownames(hess) <- names(par)
  colnames(hess) <- names(par)

  hess['c', 'c'] <- ne/c^2/2 - 1/c^3 * sum(r^2)

  hess['c', 'alp.0'] <- -1/c^2 * sum(r)
  hess['alp.0', 'c'] <- hess['c', 'alp.0']

  if(nx > 0){
    hess['c', name.alp.x] <- -1/c^2 * t(t(ex) %*% r)
    hess[name.alp.x, 'c'] <- t(hess['c', name.alp.x])
  }

  if(nh > 0){
    hess['c', name.alp.h] <- -1/c^2 * t(t(eh) %*% r)
    hess[name.alp.h, 'c'] <- t(hess['c', name.alp.h])
  }

  hess['c', name.alp.g] <- -1/c^2 * t(t(eg) %*% r)
  hess[name.alp.g, 'c'] <- t(hess['c', name.alp.g])

  if(est.b){
    hess['c', 'b'] <- -1/c^2 * sum(ed * r)
    hess['b', 'c'] <- hess['c', 'b']
  }

  ############

  hess['alp.0', 'alp.0'] <- -ne/c

  if(nx > 0){
    hess['alp.0', name.alp.x] <- -1/c * t(t(ex) %*% one)
    hess[name.alp.x, 'alp.0'] <- t(hess['alp.0', name.alp.x])
  }

  if(nh > 0){
    hess['alp.0', name.alp.h] <- -1/c * t(t(eh) %*% one)
    hess[name.alp.h, 'alp.0'] <- t(hess['alp.0', name.alp.h])
  }

  hess['alp.0', name.alp.g] <- -1/c * t(t(eg) %*% one)
  hess[name.alp.g, 'alp.0'] <- t(hess['alp.0', name.alp.g])

  if(est.b){
    hess['alp.0', 'b'] <- -1/c * sum(ed)
    hess['b', 'alp.0'] <- hess['alp.0', 'b']
  }

  ############

  if(nx > 0){
    hess[name.alp.x, name.alp.x] <- -1/c * (t(ex) %*% ex)

    if(nh > 0){
      hess[name.alp.x, name.alp.h] <- -1/c * (t(ex) %*% eh)
      hess[name.alp.h, name.alp.x] <- t(hess[name.alp.x, name.alp.h])
    }

    hess[name.alp.x, name.alp.g] <- -1/c * (t(ex) %*% eg)
    hess[name.alp.g, name.alp.x] <- t(hess[name.alp.x, name.alp.g])

    if(est.b){
      hess[name.alp.x, 'b'] <- -1/c * (t(ex) %*% ed)
      hess['b', name.alp.x] <- t(hess[name.alp.x, 'b'])
    }
  }

  ############

  if(nh > 0){
    hess[name.alp.h, name.alp.h] <- -1/c * (t(eh) %*% eh)

    hess[name.alp.h, name.alp.g] <- -1/c * (t(eh) %*% eg)
    hess[name.alp.g, name.alp.h] <- t(hess[name.alp.h, name.alp.g])

    if(est.b){
      hess[name.alp.h, 'b'] <- -1/c * (t(eh) %*% ed)
      hess['b', name.alp.h] <- t(hess[name.alp.h, 'b'])
    }
  }

  ############

  hess[name.alp.g, name.alp.g] <- -1/c * (t(eg) %*% eg)

  if(est.b){
    hess[name.alp.g, 'b'] <- -1/c * (t(eg) %*% ed)
    hess['b', name.alp.g] <- t(hess[name.alp.g, 'b'])
  }

  ############

  if(est.b){
    hess['b', 'b'] <- -1/c * sum(ed)
  }

  hess['a', 'a'] <- sum(xi)

  if(nx > 0){
    hess['a', name.bet.x] <- t(t(rx) %*% xi)
    hess[name.bet.x, 'a'] <- t(hess['a', name.bet.x])
  }

  if(ny > 0){
    hess['a', name.bet.y] <- t(t(ry) %*% xi)
    hess[name.bet.y, 'a'] <- t(hess['a', name.bet.y])
  }

  hess['a', name.bet.z] <- sum(xi * eta)
  hess[name.bet.z, 'a'] <- t(hess['a', name.bet.z])

  hess['a', name.alp.g] <- bet.z * t(t(rg) %*% xi)

  ##############

  if(nx > 0){
    hess[name.bet.x, name.bet.x] <- t(rx) %*% (xi * rx)

    if(ny > 0){
      hess[name.bet.x, name.bet.y] <- t(rx) %*% (xi * ry)
      hess[name.bet.y, name.bet.x] <- t(hess[name.bet.x, name.bet.y])
    }

    hess[name.bet.x, name.bet.z] <- t(rx) %*% (xi * eta)
    hess[name.bet.z, name.bet.x] <- t(hess[name.bet.x, name.bet.z])

    hess[name.bet.x, name.alp.g] <- bet.z * (t(rx) %*% (xi * rg))
  }

  ############

  if(ny > 0){
    hess[name.bet.y, name.bet.y] <- t(ry) %*% (xi * ry)

    hess[name.bet.y, name.bet.z] <- t(ry) %*% (xi * eta)
    hess[name.bet.z, name.bet.y] <- t(hess[name.bet.y, name.bet.z])

    hess[name.bet.y, name.alp.g] <- bet.z * (t(ry) %*% (xi * rg))
  }

  ############

  hess[name.bet.z, name.bet.z] <- sum(xi * eta^2)

  hess[name.bet.z, name.alp.g] <- t(rg) %*% (rd + Delta) + bet.z * (t(rg) %*% (xi * eta))

  hess

}

