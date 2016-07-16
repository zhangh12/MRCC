
fisher.info <- function(par, rdata, edata, par.pos){

  id3 <- intersect(rownames(rdata$data), rownames(edata$data))
  id1 <- setdiff(rownames(rdata$data), id3)
  id2 <- setdiff(rownames(edata$data), id3)

  m1 <- length(id1)
  m2 <- length(id2)
  m3 <- length(id3)

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
  names(prd) <- rownames(rdata$data)

  if(length(rdata$vy) > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
  }

  rd <- as.vector(rdata$data[, rdata$vd])
  n1 <- sum(rd)
  n0 <- m1 + m3 - n1

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
  names(Delta) <- rownames(rdata$data)

  xi <- -d.hat * (1 - d.hat)
  names(xi) <- rownames(rdata$data)

  res <- edata$data[, edata$vz] - alp.0 - eg %*% alp.g
  if(!is.na(b)){
    res <- res - (bet.z * exp(c) + b) * edata$data[, edata$vd]
  }

  if(length(edata$vx) > 0){
    res <- res - ex %*% alp.x
  }
  res <- as.vector(res)
  names(res) <- rownames(edata$data)

  ###########################
  ## calculate information ##
  ###########################

  sc <- matrix(0, nrow = length(par), ncol = m1 + m2 + m3)
  rownames(sc) <- names(par)
  colnames(sc) <- c(id1, id2, id3)

  id13 <- c(id1, id3)
  id23 <- c(id2, id3)

  ## a
  sc['a', id13] <- Delta[id13]

  ## c
  sc['c', id13] <- exp(c)/2 * bet.z^2 * Delta[id13]
  sc['c', id23] <- sc['c', id23] - 1/2 + exp(-c)/2 * res[id23]^2

  ## bet.x
  if(length(rdata$vx) > 0){
    sc[paste0('bet.', rdata$vx), id13] <- t(Delta[id13] * rx[id13, , drop = FALSE])
  }

  ## bet.y
  if(length(rdata$vy) > 0){
    sc[paste0('bet.', rdata$vy), id13] <- t(Delta[id13] * ry[id13, , drop = FALSE])
  }

  ## bet.z
  sc[paste0('bet.', edata$vz), id13] <- Delta[id13] * prd[id13]

  ## alp.0
  sc['alp.0', id13] <- bet.z * Delta[id13]
  sc['alp.0', id23] <- sc['alp.0', id23] + exp(-c) * res[id23]

  ## alp.x
  if(length(edata$vx) > 0){
    sc[paste0('alp.', rdata$vx), id13] <- bet.z * t(Delta[id13] * rx[id13, , drop = FALSE])
    sc[paste0('alp.', rdata$vx), id23] <- sc[paste0('alp.', rdata$vx), id23] + exp(-c) * t(res[id23] * ex[id23, , drop = FALSE])
  }

  ## alp.g
  sc[paste0('alp.', rdata$vg), id13] <- bet.z * t(Delta[id13] * rg[id13, , drop = FALSE])
  sc[paste0('alp.', edata$vg), id23] <- sc[paste0('alp.', edata$vg), id23] + exp(-c) * t(res[id23] * eg[id23, , drop = FALSE])

  if(length(id1) > 0){
    cov1 <- (m1-1)*cov(t(sc[, id1]))
  }else{
    cov1 <- NULL
  }

  if(length(id2) > 0){
    cov2 <- (m2-1)*cov(t(sc[, id2]))
  }else{
    cov2 <- NULL
  }

  if(length(id3) > 0){
    cov3 <- (m3-1)*cov(t(sc[, id3]))
  }else{
    cov3 <- NULL
  }

  info <- NULL

  if(!is.null(cov1)){
    if(is.null(info)){
      info <- cov1
    }else{
      info <- info + cov1
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