
# calculate hessian matrix without component of bet (bet is given)
hessian.LRT <- function(par, rdata, edata, par.pos, bet){

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

  bet.z <- bet

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

  r <- ez - alp.0 - eg %*% alp.g - (bet.z * c + b) * ed

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

  ###########################
  ## calculate Hess matrix ##
  ###########################

  hess <- matrix(0, nrow = length(par), ncol = length(par))
  rownames(hess) <- names(par)
  colnames(hess) <- names(par)

  hess['c', 'c'] <- ne/c^2/2 - 1/c^3 * sum(r^2) -
    2 * bet.z/c^2 * sum(ed * r) - bet.z^2/c * sum(ed)

  hess['c', 'alp.0'] <- -1/c^2 * sum(r) -
    bet.z/c * sum(ed)

  if(nx > 0){
    hess['c', name.alp.x] <- -1/c^2 * t(t(ex) %*% r) -
      bet.z/c * t(t(ex) %*% ed)
  }

  if(nh > 0){
    hess['c', name.alp.h] <- -1/c^2 * t(t(eh) %*% r) -
      bet.z/c * t(t(eh) %*% ed)
  }

  hess['c', name.alp.g] <- -1/c^2 * t(t(eg) %*% r) -
    bet.z/c * t(t(eg) %*% ed)

  if(est.b){
    hess['c', 'b'] <- -bet.z/c * sum(ed) - 1/c^2 * sum(ed * r)
  }

  ############

  hess['alp.0', 'alp.0'] <- -ne/c

  if(nx > 0){
    hess['alp.0', name.alp.x] <- -1/c * t(t(ex) %*% one)
  }

  if(nh > 0){
    hess['alp.0', name.alp.h] <- -1/c * t(t(eh) %*% one)
  }

  hess['alp.0', name.alp.g] <- -1/c * t(t(eg) %*% one)

  if(est.b){
    hess['alp.0', 'b'] <- -1/c * sum(ed)
  }

  ############

  if(nx > 0){
    hess[name.alp.x, name.alp.x] <- -1/c * (t(ex) %*% ex)

    if(nh > 0){
      hess[name.alp.x, name.alp.h] <- -1/c * (t(ex) %*% eh)
    }

    hess[name.alp.x, name.alp.g] <- -1/c * (t(ex) %*% eg)

    if(est.b){
      hess[name.alp.x, 'b'] <- -1/c * (t(ex) %*% ed)
    }
  }

  ############

  if(nh > 0){
    hess[name.alp.h, name.alp.h] <- -1/c * (t(eh) %*% eh)

    hess[name.alp.h, name.alp.g] <- -1/c * (t(eh) %*% eg)

    if(est.b){
      hess[name.alp.h, 'b'] <- -1/c * (t(eh) %*% ed)
    }
  }

  ############

  hess[name.alp.g, name.alp.g] <- bet.z^2 * (t(rg) %*% (xi * rg)) - 1/c * (t(eg) %*% eg)

  if(est.b){
    hess[name.alp.g, 'b'] <- -1/c * (t(eg) %*% ed)
  }

  hess[name.alp.g, 'a'] <- bet.z * (t(rg) %*% xi)

  if(nx > 0){
    hess[name.alp.g, name.bet.x] <- bet.z * (t(rg) %*% (xi * rx))
  }

  if(ny > 0){
    hess[name.alp.g, name.bet.y] <- bet.z * (t(rg) %*% (xi * ry))
  }

  ############

  if(est.b){
    hess['b', 'b'] <- -1/c * sum(ed)
  }

  hess['a', 'a'] <- sum(xi)

  if(nx > 0){
    hess['a', name.bet.x] <- t(t(rx) %*% xi)
  }

  if(ny > 0){
    hess['a', name.bet.y] <- t(t(ry) %*% xi)
  }

  if(nx > 0){
    hess[name.bet.x, name.bet.x] <- t(rx) %*% (xi * rx)

    if(ny > 0){
      hess[name.bet.x, name.bet.y] <- t(rx) %*% (xi * ry)
    }
  }

  ############

  if(ny > 0){
    hess[name.bet.y, name.bet.y] <- t(ry) %*% (xi * ry)
  }

  ############

  for(i in 1:nrow(hess)){
    for(j in 1:ncol(hess)){
      if(hess[i, j] == 0){
        hess[i, j] <- hess[j, i]
      }
    }
  }

  hess <- (hess + t(hess))/2

  #N <- length(unique(c(rownames(rdata$data)), rownames(edata$data)))
  #hess/N
  hess

}

