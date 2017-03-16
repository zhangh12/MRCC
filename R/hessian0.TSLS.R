
hessian0.TSLS <- function(par, rdata, edata, omega, par.pos){

  if(any(is.na(par))){
    return(NULL)
  }

  id1 <- omega$id1
  id10 <- omega$id10
  id11 <- omega$id11

  id2 <- omega$id2
  id20 <- omega$id20
  id21 <- omega$id21

  id3 <- omega$id3
  id30 <- omega$id30
  id31 <- omega$id31

  m1 <- length(id1)
  m10 <- length(id10)
  m11 <- length(id11)

  m2 <- length(id2)
  m20 <- length(id20)
  m21 <- length(id21)

  m3 <- length(id3)
  m30 <- length(id30)
  m31 <- length(id31)

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
  if(m10 > 0){
    rg10 <- rg[id10, , drop = FALSE]
  }
  if(m11 > 0){
    rg11 <- rg[id11, , drop = FALSE]
  }
  if(m20 > 0){
    eg20 <- eg[id20, , drop = FALSE]
  }
  if(m21 > 0){
    eg21 <- eg[id21, , drop = FALSE]
  }
  if(m30 > 0){
    rg30 <- rg[id30, , drop = FALSE]
    eg30 <- eg[id30, , drop = FALSE]
  }
  if(m31 > 0){
    rg31 <- rg[id31, , drop = FALSE]
    eg31 <- eg[id31, , drop = FALSE]
  }
  if(m2 > 0){
    eg2 <- eg[id2, , drop = FALSE]
  }
  if(m3 > 0){
    eg3 <- eg[id3, , drop = FALSE]
  }

  eta <- as.vector(rg %*% alp.g)

  if(nx > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    if(m10 > 0){
      rx10 <- rx[id10, , drop = FALSE]
    }
    if(m11 > 0){
      rx11 <- rx[id11, , drop = FALSE]
    }
    if(m20 > 0){
      ex20 <- ex[id20, , drop = FALSE]
    }
    if(m21 > 0){
      ex21 <- ex[id21, , drop = FALSE]
    }
    if(m30 > 0){
      rx30 <- rx[id30, , drop = FALSE]
      ex30 <- ex[id30, , drop = FALSE]
    }
    if(m31 > 0){
      rx31 <- rx[id31, , drop = FALSE]
      ex31 <- ex[id31, , drop = FALSE]
    }
    if(m2 > 0){
      ex2 <- ex[id2, , drop = FALSE]
    }
    if(m3 > 0){
      ex3 <- ex[id3, , drop = FALSE]
    }
  }

  if(nh > 0){
    eh <- as.matrix(edata$data[, edata$vh, drop = FALSE])

    if(m2 > 0){
      eh2 <- eh[id2, , drop = FALSE]
    }
    if(m3 > 0){
      eh3 <- eh[id3, , drop = FALSE]
    }
    if(m21 > 0){
      eh21 <- eh[id21, , drop = FALSE]
    }
    if(m31 > 0){
      eh31 <- eh[id31, , drop = FALSE]
    }
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

  p0 <- 1/ n0 / (1 + n1/n0 * delta) # 1 / (n0 + n1 * delta)
  p1 <- p0 * delta
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

  one.e <- rep(1, ne)
  one1 <- rep(1, m1)
  one2 <- rep(1, m2)
  one3 <- rep(1, m3)
  one21 <- rep(1, m21)
  one31 <- rep(1, m31)
  one.r <- rep(1, nr)

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

  hess['c', 'c'] <- -ne/c^2/2

  ############

  hess['alp.0', 'alp.0'] <- -ne/c

  if(nx > 0){
    if(m30 > 0){
      hess['alp.0', name.alp.x] <- -1/c * (m30 * t(t(rx) %*% p0))
    }

    if(m31 > 0){
      hess['alp.0', name.alp.x] <- hess['alp.0', name.alp.x] - 1/c * (m31 * t(t(rx) %*% p1))
    }

    if(m2 > 0){
      hess['alp.0', name.alp.x] <- hess['alp.0', name.alp.x] - 1/c * (t(ex2) %*% one2)
    }

    hess[name.alp.x, 'alp.0'] <- t(hess['alp.0', name.alp.x])
  }

  if(nh > 0){
    hess['alp.0', name.alp.h] <- -1/c * t(t(eh) %*% one.e)

    hess[name.alp.h, 'alp.0'] <- t(hess['alp.0', name.alp.h])
  }

  if(m30 > 0){
    hess['alp.0', name.alp.g] <- -1/c * (m30 * (t(rg) %*% p0))
  }
  if(m31 > 0){
    hess['alp.0', name.alp.g] <- hess['alp.0', name.alp.g] - 1/c * (m31 * (t(rg) %*% p1))
  }
  if(m2 > 0){
    hess['alp.0', name.alp.g] <- hess['alp.0', name.alp.g] - 1/c * (t(eg2) %*% one2)
  }
  hess[name.alp.g, 'alp.0'] <- t(hess['alp.0', name.alp.g])

  if(est.b){
    hess['alp.0', 'b'] <- -1/c * sum(ed)
    hess['b', 'alp.0'] <- hess['alp.0', 'b']
  }

  ############

  if(nx > 0){
    if(m30 > 0){
      hess[name.alp.x, name.alp.x] <- -1/c * (m30 * (t(p0 * rx) %*% rx))
    }
    if(m31 > 0){
      hess[name.alp.x, name.alp.x] <- hess[name.alp.x, name.alp.x] - 1/c * (m31 * (t(p1 * rx) %*% rx))
    }
    if(m2 > 0){
      hess[name.alp.x, name.alp.x] <- hess[name.alp.x, name.alp.x] - 1/c * (t(ex2) %*% ex2)
    }

    if(nh > 0){
      hess[name.alp.x, name.alp.h] <- -1/c * (t(ex) %*% eh)
      hess[name.alp.h, name.alp.x] <- t(hess[name.alp.x, name.alp.h])
    }

    if(m30 > 0){
      hess[name.alp.x, name.alp.g] <- -1/c * (m30 * (t(p0 * rx) %*% rg))
    }
    if(m31 > 0){
      hess[name.alp.x, name.alp.g] <- hess[name.alp.x, name.alp.g] - 1/c * (m31 * (t(p1 * rx) %*% rg))
    }
    if(m2 > 0){
      hess[name.alp.x, name.alp.g] <- hess[name.alp.x, name.alp.g] - 1/c * (t(ex2) %*% eg2)
    }
    hess[name.alp.g, name.alp.x] <- t(hess[name.alp.x, name.alp.g])

    if(est.b){
      if(m31 > 0){
        hess[name.alp.x, 'b'] <- -1/c * (m31 * t(t(rx) %*% p1))
      }
      if(m21 > 0){
        hess[name.alp.x, 'b'] <- hess[name.alp.x, 'b'] - 1/c * (t(ex21) %*% one21)
      }
      hess['b', name.alp.x] <- t(hess[name.alp.x, 'b'])
    }
  }

  ############

  if(nh > 0){
    hess[name.alp.h, name.alp.h] <- -1/c * (t(eh) %*% eh)

    hess[name.alp.h, name.alp.g] <- -1/c * (t(eh) %*% eg)

    hess[name.alp.g, name.alp.h] <- t(hess[name.alp.h, name.alp.g])

    if(est.b){
      if(m31 > 0){
        hess[name.alp.h, 'b'] <- -1/c * (t(eh31) %*% one31)
      }
      if(m21 > 0){
        hess[name.alp.h, 'b'] <- hess[name.alp.h, 'b'] - 1/c * (t(eh21) %*% one21)
      }
      hess['b', name.alp.h] <- t(hess[name.alp.h, 'b'])
    }
  }

  ############

  if(m30 > 0){
    hess[name.alp.g, name.alp.g] <- - 1/c * (m30 * (t(p0 * rg) %*% rg))
  }
  if(m31 > 0){
    hess[name.alp.g, name.alp.g] <- hess[name.alp.g, name.alp.g] - 1/c * (m31 * (t(p1 * rg) %*% rg))
  }
  if(m2 > 0){
    hess[name.alp.g, name.alp.g] <- hess[name.alp.g, name.alp.g] - 1/c * (t(eg2) %*% eg2)
  }

  if(est.b){
    if(m31 > 0){
      hess[name.alp.g, 'b'] <- -1/c * (m31 * (t(rg) %*% p1))
    }
    if(m21 > 0){
      hess[name.alp.g, 'b'] <- hess[name.alp.g, 'b'] - 1/c * (t(eg21) %*% one21)
    }
    hess['b', name.alp.g] <- t(hess[name.alp.g, 'b'])
  }

  ############

  if(est.b){
    hess['b', 'b'] <- -1/c * sum(ed)
  }

  ############

  hess['a', 'a'] <- n0 * sum(p0 * Delta)

  if(nx > 0){
    hess['a', name.bet.x] <- n0 * t(t(rx) %*% (p0 * Delta))
    hess[name.bet.x, 'a'] <- t(hess['a', name.bet.x])
  }

  if(ny > 0){
    hess['a', name.bet.y] <- n0 * t(t(ry) %*% (p0 * Delta))
    hess[name.bet.y, 'a'] <- t(hess['a', name.bet.y])
  }

  hess['a', name.bet.z] <- n0 * sum(p0 * Delta * eta)
  hess[name.bet.z, 'a'] <- hess['a', name.bet.z]

  hess['a', name.alp.g] <- bet.z * n0 * t(t(rg) %*% (p0 * Delta))

  ##########

  if(nx > 0){
    hess[name.bet.x, name.bet.x] <- n0 * (t(rx) %*% (p0 * Delta * rx))

    if(ny > 0){
      hess[name.bet.x, name.bet.y] <- n0 * (t(rx) %*% (p0 * Delta * ry))
      hess[name.bet.y, name.bet.x] <- t(hess[name.bet.x, name.bet.y])
    }

    hess[name.bet.x, name.bet.z] <- n0 * (t(rx) %*% (p0 * Delta * eta))
    hess[name.bet.z, name.bet.x] <- t(hess[name.bet.x, name.bet.z])

    hess[name.bet.x, name.alp.g] <- bet.z * n0 * (t(rx) %*% (p0 * Delta * rg))
  }

  ############

  if(ny > 0){
    hess[name.bet.y, name.bet.y] <- n0 * (t(ry) %*% (p0 * Delta * ry))

    hess[name.bet.y, name.bet.z] <- n0 * (t(ry) %*% (p0 * Delta * eta))
    hess[name.bet.z, name.bet.y] <- t(hess[name.bet.y, name.bet.z])

    hess[name.bet.y, name.alp.g] <- bet.z * n0 * (t(ry) %*% (p0 * Delta * rg))
  }

  ############

  hess[name.bet.z, name.bet.z] <- n0 * sum(p0 * Delta * eta^2)

  hess[name.bet.z, name.alp.g] <- bet.z * n0 * t(t(rg) %*% (p0 * Delta * eta))

  hess

}

