

score <- function(par, data.retro, data.expo, par.pos){

  ne <- nrow(data.expo$edat)
  nr <- nrow(data.retro$rdat)

  a <- par['a']

  if('b' %in% par.pos$par.name){
    b <- par['b']
  }else{
    b <- NA
  }

  c <- par['c']
  alp.0 <- par['alp.0']

  if('alp.o' %in% par.pos$par.name){
    alp.o <- par[par.pos['alp.o', 'start']:par.pos['alp.o', 'end']]
  }else{
    alp.o <- NA
  }

  alp.g <- par[par.pos['alp.g', 'start']:par.pos['alp.g', 'end']]

  if('bet.a' %in% par.pos$par.name){
    bet.a <- par[par.pos['bet.a', 'start']:par.pos['bet.a', 'end']]
  }else{
    bet.a <- NA
  }

  if('bet.o' %in% par.pos$par.name){
    bet.o <- par[par.pos['bet.o', 'start']:par.pos['bet.o', 'end']]
  }else{
    bet.o <- NA
  }

  bet.e <- par[par.pos['bet.e', 'start']:par.pos['bet.e', 'end']]

  yr <- as.vector(data.retro$rdat[, data.retro$retro.var])
  n1 <- sum(yr)
  n0 <- nr - n1

  const <- a + alp.0 * bet.e + .5 * exp(c) * bet.e^2
  if(!is.na(b)){
    const <- const + b * bet.e
  }

  lin <- as.matrix(data.retro$rdat[, data.retro$geno.var, drop = FALSE]) %*% alp.g * bet.e
  if(length(data.retro$overlap.covar) > 0){
    lin <- lin + as.matrix(data.retro$rdat[, data.retro$overlap.covar, drop = FALSE]) %*% (bet.o + bet.e * alp.o)
  }

  if(length(data.retro$add.covar.retro) > 0){
    lin <- lin + as.matrix(data.retro$rdat[, data.retro$add.covar.retro, drop = FALSE]) %*% bet.a
  }

  lin <- const + lin
  lin <- as.vector(lin)

  delta <- exp(lin)
  p1 <- 1 - 1/(1 + n1 / n0 * delta)
  d <- p1 * (1-p1)

  res <- data.expo$edat[, data.expo$expo.var] - alp.0 - as.matrix(data.expo$edat[, data.expo$geno.var]) %*% alp.g
  if(!is.na(b)){
    res <- res - (bet.e * exp(c) + b) * data.expo$edat[, data.expo$retro.var]
  }

  if(length(data.expo$overlap.covar) > 0){
    res <- res - as.matrix(data.expo$edat[, data.expo$overlap.covar, drop = FALSE]) %*% alp.o
  }

  var.hat <- mean(res^2)

  rg <- as.matrix(data.retro$rdat[, data.retro$geno.var, drop = FALSE])
  eg <- as.matrix(data.expo$edat[, data.expo$geno.var, drop = FALSE])

  if(length(data.retro$overlap.covar) > 0){
    ro <- as.matrix(data.retro$rdat[, data.retro$overlap.covar, drop = FALSE])
    eo <- as.matrix(data.expo$edat[, data.expo$overlap.covar, drop = FALSE])
  }

  if(length(data.retro$add.covar.retro) > 0){
    ra <- as.matrix(data.retro$rdat[, data.retro$add.covar.retro, drop = FALSE])
  }

  ############################
  ## calculate score vector ##
  ############################

  gr <- rep(NA, length(par))
  names(gr) <- names(par)

  ell.a <- sum(yr - p1)
  gr['a'] <- ell.a

  ell.c <- .5 * exp(c) * bet.e^2 * ell.a - .5 * ne * (1 - exp(-c) * var.hat)
  gr['c'] <- ell.c

  ell.alp.0 <- bet.e * ell.a + exp(-c) * sum(res)
  gr['alp.0'] <- ell.alp.0

  if(length(data.retro$overlap.covar) > 0){
    ell.alp.o <- bet.e * t(ro) %*% (yr - p1) + exp(-c) * t(eo) %*% res
    ell.alp.o <- as.vector(ell.alp.o)
    gr[paste0('alp.', data.expo$overlap.covar)] <- ell.alp.o

    ell.bet.o <- t(ro) %*% (yr - p1)
    ell.bet.o <- as.vector(ell.bet.o)
    gr[paste0('bet.', data.retro$overlap.covar)] <- ell.bet.o
  }

  ell.alp.g <- bet.e * t(rg) %*% (yr - p1) + exp(-c) * t(eg) %*% res
  ell.alp.g <- as.vector(ell.alp.g)
  gr[paste0('alp.', data.expo$geno.var)] <- ell.alp.g

  if(length(data.retro$add.covar.retro) > 0){
    ell.bet.a <- t(ra) %*% (yr - p1)
    ell.bet.a <- as.vector(ell.bet.a)
    gr[paste0('bet.', data.retro$add.covar.retro)] <- ell.bet.a
  }

  ell.bet.e <- (alp.0 + exp(c) * bet.e) * ell.a + sum((yr - p1) * (rg %*% alp.g))
  if(length(data.retro$overlap.covar) > 0){
    ell.bet.e <- ell.bet.e + sum((yr - p1) * (ro %*% alp.o))
  }
  gr[paste0('bet.', data.expo$expo.var)] <- ell.bet.e

  gr <- - gr

  gr


}