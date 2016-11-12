
fisher.info0.TSLS <- function(par, rdata, edata, par.pos){

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

  n0 <- n10 + m3
  n1 <- n11

  ########
  if('alp.x' %in% par.pos$par.name){
    alp.x <- par[par.pos['alp.x', 'start']:par.pos['alp.x', 'end']]
    nx <- length(alp.x)
    rm(alp.x)
  }else{
    nx <- 0
  }

  alp.g <- par[par.pos['alp.g', 'start']:par.pos['alp.g', 'end']]
  ng <- length(alp.g)
  rm(alp.g)

  if('bet.y' %in% par.pos$par.name){
    bet.y <- par[par.pos['bet.y', 'start']:par.pos['bet.y', 'end']]
    ny <- length(bet.y)
    rm(bet.y)
  }else{
    ny <- 0
  }

  name <- c('a')
  name.alp.g <- paste0('alp.', rdata$vg)
  if(nx > 0){
    name.alp.x <- paste0('alp.', rdata$vx)
    name.bet.x <- paste0('bet.', rdata$vx)
    name <- c(name, name.bet.x)
  }

  if(ny > 0){
    name.bet.y <- paste0('bet.', rdata$vy)
    name <- c(name, name.bet.y)
  }

  name.bet.z <- paste0('bet.', edata$vz)
  name <- c(name, name.bet.z)

  info <- -hessian0.TSLS(par, rdata, edata, par.pos)

  info[name, name.alp.g] <- 0
  C1 <- info[name, 'a', drop = FALSE]
  info[name, name] <- info[name, name] - (n0+n1)/n0/n1 * (C1 %*% t(C1))
  info


}
