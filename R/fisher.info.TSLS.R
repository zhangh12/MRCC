
fisher.info.TSLS <- function(par, rdata, edata, omega, par.pos){

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
  names(eta) <- rownames(rdata$data)

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
  names(rd) <- rownames(rdata$data)
  n1 <- sum(rd)
  n0 <- nr - n1

  ed <- as.vector(edata$data[, edata$vd])
  names(ed) <- rownames(edata$data)
  ez <- as.vector(edata$data[, edata$vz])
  names(ez) <- rownames(edata$data)

  lin <- a + rg %*% alp.g * bet.z
  if(nx > 0){
    lin <- lin + rx %*% bet.x
  }

  if(ny > 0){
    lin <- lin + ry %*% bet.y
  }

  lin <- as.vector(lin)
  names(lin) <- rownames(rdata$data)

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
  names(r) <- rownames(edata$data)

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
  ## calculate information ##
  ###########################

  sc <- matrix(0, nrow = m1 + m2 + m3, ncol = length(par))
  rownames(sc) <- c(id1, id2, id3)
  colnames(sc) <- names(par)

  id.r <- c(id1, id3)
  id.e <- c(id2, id3)

  ## c
  sc[id.e, 'c'] <- -1/c/2 + 1/c^2/2 * r[id.e]^2

  ## alp.0
  sc[id.e, 'alp.0'] <- 1/c * r[id.e]

  ## alp.x
  if(nx > 0){
    sc[id.e, name.alp.x] <- 1/c * (r[id.e] * ex[id.e, , drop = FALSE])
  }

  ## alp.h
  if(nh > 0){
    sc[id.e, name.alp.h] <- 1/c * (r[id.e] * eh[id.e, , drop = FALSE])
  }

  ## alp.g
  sc[id.e, name.alp.g] <- 1/c * (r[id.e] * eg[id.e, , drop = FALSE])

  ##
  if(est.b){
    sc[id.e, 'b'] <- 1/c * ed[id.e] * r[id.e]
  }

  ## a
  sc[id.r, 'a'] <- rd[id.r] + Delta[id.r]

  ## bet.x
  if(nx > 0){
    sc[id.r, name.bet.x] <- (rd[id.r] + Delta[id.r]) * rx[id.r, , drop = FALSE]
  }

  ## bet.y
  if(ny > 0){
    sc[id.r, name.bet.y] <- (rd[id.r] + Delta[id.r]) * ry[id.r, , drop = FALSE]
  }

  ## bet.z
  sc[id.r, name.bet.z] <- (rd[id.r] + Delta[id.r]) * eta[id.r]

  n <- m1 + m2 + m3

  if(m10 > 0){
    cov10 <- (m10 - 1) * cov(sc[id10, ])
  }else{
    cov10 <- NULL
  }

  if(m11 > 0){
    cov11 <- (m11 - 1) * cov(sc[id11, ])
  }else{
    cov11 <- NULL
  }

  if(m20 > 0){
    cov20 <- (m20 - 1) * cov(sc[id20, ])
  }else{
    cov20 <- NULL
  }

  if(m21 > 0){
    cov21 <- (m21 - 1) * cov(sc[id21, ])
  }else{
    cov21 <- NULL
  }

  if(m30 > 0){
    cov30 <- (m30 - 1) * cov(sc[id30, ])
  }else{
    cov30 <- NULL
  }

  if(m31 > 0){
    cov31 <- (m31 - 1) * cov(sc[id31, ])
  }else{
    cov31 <- NULL
  }

  info <- NULL

  if(!is.null(cov10)){
    if(is.null(info)){
      info <- cov10
    }else{
      info <- info + cov10
    }
  }

  if(!is.null(cov11)){
    if(is.null(info)){
      info <- cov11
    }else{
      info <- info + cov11
    }
  }

  if(!is.null(cov20)){
    if(is.null(info)){
      info <- cov20
    }else{
      info <- info + cov20
    }
  }

  if(!is.null(cov21)){
    if(is.null(info)){
      info <- cov21
    }else{
      info <- info + cov21
    }
  }

  if(!is.null(cov30)){
    if(is.null(info)){
      info <- cov30
    }else{
      info <- info + cov30
    }
  }

  if(!is.null(cov31)){
    if(is.null(info)){
      info <- cov31
    }else{
      info <- info + cov31
    }
  }

  info

}