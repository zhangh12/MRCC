
hessian0.TSLS <- function(par, rdata, edata, par.pos){

  if(any(is.na(par))){
    return(NULL)
  }

  id3 <- intersect(rownames(rdata$data), rownames(edata$data))
  id1 <- setdiff(rownames(rdata$data), id3)
  id2 <- setdiff(rownames(edata$data), id3)

  id10 <- setdiff(rownames(rdata$data[rdata$data[, rdata$vd] == 0, , drop = FALSE]), id3)
  id11 <- setdiff(rownames(rdata$data[rdata$data[, rdata$vd] == 1, , drop = FALSE]), id3)

  m1 <- length(id1)
  m2 <- length(id2)
  m3 <- length(id3)

  n10 <- length(id10)
  n11 <- length(id11)

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
  eg2 <- eg[id2, , drop = FALSE]
  eta <- as.vector(rg %*% alp.g)

  if(nx > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    ex2 <- ex[id2, , drop = FALSE]
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
  one2 <- rep(1, m2)
  one13 <- rep(1, nr)

  name.alp.g <- paste0('alp.', rdata$vg)
  if(nx > 0){
    name.alp.x <- paste0('alp.', rdata$vx)
    name.bet.x <- paste0('bet.', rdata$vx)
  }

  if(ny > 0){
    name.bet.y <- paste0('bet.', rdata$vy)
  }

  name.bet.z <- paste0('bet.', edata$vz)


  ###########################
  ## calculate Hess matrix ##
  ###########################

  hess <- matrix(0, nrow = length(par), ncol = length(par))
  rownames(hess) <- names(par)
  colnames(hess) <- names(par)

  hess['c', 'c'] <- -ne/c^2/2

  ############

  hess['alp.0', 'alp.0'] <- -ne/c

  if(nx > 0){
    hess['alp.0', name.alp.x] <- -1/c * (t(t(ex2) %*% one2) + m3 * t(t(p * rx) %*% one13))
    hess[name.alp.x, 'alp.0'] <- t(hess['alp.0', name.alp.x, drop = FALSE])
  }

  hess['alp.0', name.alp.g] <- -1/c * (t(t(eg2) %*% one2) + m3 * t(t(p * rg) %*% one13))
  hess[name.alp.g, 'alp.0'] <- t(hess['alp.0', name.alp.g, drop = FALSE])

  ############

  if(nx > 0){
    hess[name.alp.x, name.alp.x] <- -1/c * (t(t(ex2) %*% ex2) + m3 * t(t(p * rx) %*% rx))

    hess[name.alp.x, name.alp.g] <- -1/c * (t(t(ex2) %*% eg2) + m3 * t(t(p * rx) %*% rg))
    hess[name.alp.g, name.alp.x] <- t(hess[name.alp.x, name.alp.g, drop = FALSE])
  }

  ############

  #a33 <- -bet.z^2 * n0 * (t(rg) %*% (p * Delta * rg))
  hess[name.alp.g, name.alp.g] <- - 1/c * (t(t(eg2) %*% eg2) + m3 * t(t(p * rg) %*% rg))

  hess['a', name.alp.g] <- bet.z * n0 * t(t(rg) %*% (p * Delta))

  if(nx > 0){
    hess[name.bet.x, name.alp.g] <- bet.z * n0 * t(t(rg) %*% (p * Delta * rx))
  }

  if(ny > 0){
    hess[name.bet.y, name.alp.g] <- bet.z * n0 * t(t(rg) %*% (p * Delta * ry))
  }

  hess[name.bet.z, name.alp.g] <- bet.z * n0 * t(t(rg) %*% (p * Delta * eta))

  ############

  hess['a', 'a'] <- n0 * sum(p * Delta)

  if(nx > 0){
    hess['a', name.bet.x] <- n0 * t(t(rx) %*% (p * Delta))
    hess[name.bet.x, 'a'] <- t(hess['a', name.bet.x, drop = FALSE])
  }

  if(ny > 0){
    hess['a', name.bet.y] <- n0 * t(t(ry) %*% (p * Delta))
    hess[name.bet.y, 'a'] <- t(hess['a', name.bet.y, drop = FALSE])
  }

  hess['a', name.bet.z] <- n0 * sum(p * Delta * eta)
  hess[name.bet.z, 'a'] <- t(hess['a', name.bet.z, drop = FALSE])

  ##########

  if(nx > 0){
    hess[name.bet.x, name.bet.x] <- n0 * (t(rx) %*% (p * Delta * rx))

    if(ny > 0){
      hess[name.bet.x, name.bet.y] <- n0 * (t(rx) %*% (p * Delta * ry))
      hess[name.bet.y, name.bet.x] <- t(hess[name.bet.x, name.bet.y, drop = FALSE])
    }

    hess[name.bet.x, name.bet.z] <- n0 * (t(rx) %*% (p * Delta * eta))
    hess[name.bet.z, name.bet.x] <- t(hess[name.bet.x, name.bet.z, drop = FALSE])
  }

  ############

  if(ny > 0){
    hess[name.bet.y, name.bet.y] <- n0 * (t(ry) %*% (p * Delta * ry))

    hess[name.bet.y, name.bet.z] <- n0 * (t(ry) %*% (p * Delta * eta))
    hess[name.bet.z, name.bet.y] <- t(hess[name.bet.y, name.bet.z, drop = FALSE])
  }

  ############

  hess[name.bet.z, name.bet.z] <- n0 * sum(p * Delta * eta^2)

  hess

}

