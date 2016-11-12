
fisher.info <- function(par, rdata, edata, par.pos){

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
  eta <- as.vector(rg %*% alp.g)
  names(eta) <- rownames(rg)

  if(nx > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    #eta.x <- as.vector(rx %*% alp.x)
    #names(eta.x) <- rownames(rx)
  }

  if(ny > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
  }

  rd <- as.vector(rdata$data[, rdata$vd])
  names(rd) <- rownames(rdata$data)
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
  names(lin) <- rownames(rg)

  delta <- exp(lin)

  p <- 1/ n0 / (1 + n1/n0 * delta) # 1 / (n0 + n1 * delta)
  Delta <- -1 + 1 / (1 + n1/n0 * delta) # -n1 p delta
  xi <- Delta * (1 + Delta)

  r <- edata$data[, edata$vz] - alp.0 - eg %*% alp.g

  if(nx > 0){
    r <- r - ex %*% alp.x
  }
  r <- as.vector(r)
  names(r) <- rownames(eg)

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
  ## calculate information ##
  ###########################

  sc <- matrix(0, nrow = length(par), ncol = m1 + m2 + m3)
  rownames(sc) <- names(par)
  colnames(sc) <- c(id1, id2, id3)

  id13 <- c(id1, id3)
  id23 <- c(id2, id3)

  ## c
  sc['c', id23] <- -1/c/2 + 1/c^2/2 * r[id23]^2

  ## alp.0
  sc['alp.0', id23] <- 1/c * r[id23]

  ## alp.x
  if(nx > 0){
    sc[name.alp.x, id23] <- 1/c * t(r[id23] * ex[id23, , drop = FALSE])
  }

  ## alp.g
  sc[name.alp.g, id13] <- bet.z * t((rd[id13] + Delta[id13]) * rg[id13, , drop = FALSE])
  sc[name.alp.g, id23] <- sc[name.alp.g, id23] + 1/c * t(r[id23] * eg[id23, , drop = FALSE])

  ## a
  sc['a', id13] <- rd[id13] + Delta[id13]

  ## bet.x
  if(nx > 0){
    sc[name.bet.x, id13] <- t((rd[id13] + Delta[id13]) * rx[id13, , drop = FALSE])
  }

  ## bet.y
  if(ny > 0){
    sc[name.bet.y, id13] <- t((rd[id13] + Delta[id13]) * ry[id13, , drop = FALSE])
  }

  ## bet.z
  sc[name.bet.z, id13] <- (rd[id13] + Delta[id13]) * eta[id13]

  n <- length(c(id1,id2,id3))

  if(n10 > 0){
    cov10 <- (n10 - 1) * cov(t(sc[, id10]))
  }else{
    cov10 <- NULL
  }

  if(n11 > 0){
    cov11 <- (n11 - 1) * cov(t(sc[, id11]))
  }else{
    cov11 <- NULL
  }

  if(m2 > 0){
    cov2 <- (m2 - 1) * cov(t(sc[, id2]))
  }else{
    cov2 <- NULL
  }

  if(m3 > 0){
    cov3 <- (m3 - 1) * cov(t(sc[, id3]))
  }else{
    cov3 <- NULL
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

  if(!is.null(cov2)){
    if(is.null(info)){
      info <- cov2
    }else{
      info <- info + cov2
    }
  }

  if(!is.null(cov3)){
    if(is.null(info)){
      info <- cov3
    }else{
      info <- info + cov3
    }
  }

  info

}