

score.hess <- function(par, data.retro, data.expo, par.pos){

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

  ##################

  rg <- as.matrix(data.retro$rdat[, data.retro$geno.var, drop = FALSE])
  eg <- as.matrix(data.expo$edat[, data.expo$geno.var, drop = FALSE])
  prd.rg <- as.vector(rg %*% alp.g)

  if(length(data.retro$overlap.covar) > 0){
    ro <- as.matrix(data.retro$rdat[, data.retro$overlap.covar, drop = FALSE])
    eo <- as.matrix(data.expo$edat[, data.expo$overlap.covar, drop = FALSE])
    prd.ro <- as.vector(ro %*% alp.o)
  }

  if(length(data.retro$add.covar.retro) > 0){
    ra <- as.matrix(data.retro$rdat[, data.retro$add.covar.retro, drop = FALSE])
  }

  yr <- as.vector(data.retro$rdat[, data.retro$retro.var])
  n1 <- sum(yr)
  n0 <- nr - n1

  const <- a + alp.0 * bet.e + .5 * exp(c) * bet.e^2
  if(!is.na(b)){
    const <- const + b * bet.e
  }

  lin <- rg %*% alp.g * bet.e
  if(length(data.retro$overlap.covar) > 0){
    lin <- lin + ro %*% (bet.o + bet.e * alp.o)
  }

  if(length(data.retro$add.covar.retro) > 0){
    lin <- lin + ra %*% bet.a
  }

  lin <- const + lin
  lin <- as.vector(lin)

  delta <- exp(lin)
  p1 <- 1 - 1/(1 + n1 / n0 * delta)
  d <- p1 * (1-p1)

  res <- data.expo$edat[, data.expo$expo.var] - alp.0 - eg %*% alp.g
  if(!is.na(b)){
    res <- res - (bet.e * exp(c) + b) * data.expo$edat[, data.expo$retro.var]
  }

  if(length(data.expo$overlap.covar) > 0){
    res <- res - eo %*% alp.o
  }

  var.hat <- mean(res^2)


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

  ell.alp.g <- bet.e * t(rg) %*% (yr - p1) + exp(-c) * t(eg) %*% res
  ell.alp.g <- as.vector(ell.alp.g)
  gr[paste0('alp.', data.expo$geno.var)] <- ell.alp.g

  ell.bet.e <- (alp.0 + exp(c) * bet.e) * ell.a + sum((yr - p1) * prd.rg)
  gr[paste0('bet.', data.expo$expo.var)] <- ell.bet.e

  if(length(data.retro$overlap.covar) > 0){
    ell.alp.o <- bet.e * t(ro) %*% (yr - p1) + exp(-c) * t(eo) %*% res
    ell.alp.o <- as.vector(ell.alp.o)
    gr[paste0('alp.', data.expo$overlap.covar)] <- ell.alp.o

    ell.bet.o <- t(ro) %*% (yr - p1)
    ell.bet.o <- as.vector(ell.bet.o)
    gr[paste0('bet.', data.retro$overlap.covar)] <- ell.bet.o

    ell.bet.e <- ell.bet.e + sum((yr - p1) * prd.ro)
    gr[paste0('bet.', data.expo$expo.var)] <- ell.bet.e
  }

  if(length(data.retro$add.covar.retro) > 0){
    ell.bet.a <- t(ra) %*% (yr - p1)
    ell.bet.a <- as.vector(ell.bet.a)
    gr[paste0('bet.', data.retro$add.covar.retro)] <- ell.bet.a
  }

  gr <- - gr

  ###########################
  ## calculate Hess matrix ##
  ###########################

  hess <- matrix(NA, nrow = length(par), ncol = length(par))
  rownames(hess) <- names(par)
  colnames(hess) <- names(par)

  # can be defined anyway
  ell.a.a <- -sum(d)
  hess['a', 'a'] <- ell.a.a

  ell.a.c <- .5 * exp(c) * bet.e^2 * ell.a.a
  hess['a', 'c'] <- ell.a.c

  ell.a.alp.0 <- bet.e * ell.a.a
  hess['a', 'alp.0'] <- ell.a.alp.0

  ell.a.alp.g <- - bet.e * t(rg) %*% d
  ell.a.alp.g <- matrix(as.vector(ell.a.alp.g), nrow = 1)
  hess['a', paste0('alp.', data.expo$geno.var)] <- ell.a.alp.g

  ell.c.c <- matrix(.5 * exp(c) * bet.e^2 * ell.a + .5 * exp(c) * bet.e^2 * ell.a.c - .5 * ne * exp(-c) * var.hat)
  hess['c', 'c'] <- ell.c.c

  ell.c.alp.0 <- matrix(.5 * exp(c) * bet.e^2 * ell.a.alp.0 - exp(-c) * sum(res))
  hess['c', 'alp.0'] <- ell.c.alp.0

  ell.c.alp.g <- .5 * exp(c) * bet.e^2 * ell.a.alp.g - exp(-c) * t(t(eg) %*% res)
  ell.c.alp.g <- matrix(as.vector(ell.c.alp.g), nrow = 1)
  hess['c', paste0('alp.', data.expo$geno.var)] <- ell.c.alp.g

  ell.alp.0.alp.0 <- bet.e * ell.a.alp.0 - ne * exp(-c)
  hess['alp.0', 'alp.0'] <- ell.alp.0.alp.0

  ell.alp.0.alp.g <- bet.e * ell.a.alp.g - exp(-c) * colSums(eg)
  ell.alp.0.alp.g <- matrix(as.vector(ell.alp.0.alp.g), nrow = 1)
  hess['alp.0', paste0('alp.', data.expo$geno.var)] <- ell.alp.0.alp.g

  ell.alp.g.alp.g <- - t(rg) %*% (d * rg) * bet.e^2 - exp(-c) * t(eg) %*% eg
  ell.alp.g.alp.g <- as.matrix(ell.alp.g.alp.g)
  hess[paste0('alp.', data.expo$geno.var), paste0('alp.', data.expo$geno.var)] <- ell.alp.g.alp.g

  # partially depends on covariate
  ell.a.bet.e <- (alp.0 + exp(c) * bet.e) * ell.a.a - sum(prd.rg * d)
  ell.a.bet.e <- matrix(ell.a.bet.e)
  hess['a', paste0('bet.', data.expo$expo.var)] <- ell.a.bet.e

  ell.alp.g.bet.e <- t(rg) %*% (yr - p1) - t(rg) %*% d * (alp.0 + exp(c) * bet.e) * bet.e - t(rg) %*% (prd.rg * d) * bet.e
  ell.alp.g.bet.e <- matrix(as.vector(ell.alp.g.bet.e), ncol = 1)
  hess[paste0('alp.', data.expo$geno.var), paste0('bet.', data.expo$expo.var)] <- ell.alp.g.bet.e

  ell.bet.e.bet.e <- exp(c) * ell.a - sum(prd.rg * d) * (alp.0 + exp(c) * bet.e) - sum(d * prd.rg^2)
  ell.bet.e.bet.e <- matrix(ell.bet.e.bet.e)
  hess[paste0('bet.', data.expo$expo.var), paste0('bet.', data.expo$expo.var)] <- ell.bet.e.bet.e

  if(length(data.retro$overlap.covar) > 0){
    ell.a.bet.o <- - t(ro) %*% d
    ell.a.bet.o <- matrix(as.vector(ell.a.bet.o), nrow = 1)
    hess['a', paste0('bet.', data.retro$overlap.covar)] <- ell.a.bet.o

    ell.a.alp.o <- bet.e * ell.a.bet.o
    ell.a.alp.o <- matrix(as.vector(ell.a.alp.o), nrow = 1)
    hess['a', paste0('alp.', data.expo$overlap.covar)] <- ell.a.alp.o

    ell.a.bet.e <- ell.a.bet.e - sum(prd.ro * d)
    ell.a.bet.e <- matrix(ell.a.bet.e)
    hess['a', paste0('bet.', data.expo$expo.var)] <- ell.a.bet.e

    ell.c.alp.o <- .5 * exp(c) * bet.e^2 * ell.a.alp.o - exp(-c) * t(t(eo) %*% res)
    ell.c.alp.o <- matrix(as.vector(ell.c.alp.o), nrow = 1)
    hess['c', paste0('alp.', data.expo$overlap.covar)] <- ell.c.alp.o

    ell.c.bet.o <- .5 * exp(c) * bet.e^2 * ell.a.bet.o
    ell.c.bet.o <- matrix(as.vector(ell.c.bet.o), nrow = 1)
    hess['c', paste0('bet.', data.retro$overlap.covar)] <- ell.c.bet.o

    ell.alp.0.alp.o <- bet.e * ell.a.alp.o - exp(-c) * colSums(eo)
    ell.alp.0.alp.o <- matrix(as.vector(ell.alp.0.alp.o), nrow = 1)
    hess['alp.0', paste0('alp.', data.expo$overlap.covar)] <- ell.alp.0.alp.o

    ell.alp.0.bet.o <- bet.e * ell.a.bet.o
    hess['alp.0', paste0('bet.', data.retro$overlap.covar)] <- ell.alp.0.bet.o

    ell.alp.o.alp.o <- - t(ro) %*% (d * ro) * bet.e^2 - exp(-c) * t(eo) %*% eo
    ell.alp.o.alp.o <- as.matrix(ell.alp.o.alp.o)
    hess[paste0('alp.', data.expo$overlap.covar), paste0('alp.', data.expo$overlap.covar)] <- ell.alp.o.alp.o

    ell.alp.o.alp.g <- - t(ro) %*% (d * rg) * bet.e^2 - exp(-c) * t(eo) %*% eg
    ell.alp.o.alp.g <- as.matrix(ell.alp.o.alp.g)
    hess[paste0('alp.', data.expo$overlap.covar), paste0('alp.', data.expo$geno.var)] <- ell.alp.o.alp.g

    ell.alp.o.bet.o <- - t(ro) %*% (d * ro) * bet.e
    ell.alp.o.bet.o <- as.matrix(ell.alp.o.bet.o)
    hess[paste0('alp.', data.expo$overlap.covar), paste0('bet.', data.retro$overlap.covar)] <- ell.alp.o.bet.o

    ell.alp.o.bet.e <- t(ro) %*% (yr - p1) - t(ro) %*% d * (alp.0 + exp(c) * bet.e) * bet.e - t(ro) %*% (d * (prd.ro + prd.rg)) * bet.e
    ell.alp.o.bet.e <- matrix(as.vector(ell.alp.o.bet.e), ncol = 1)
    hess[paste0('alp.', data.expo$overlap.covar), paste0('bet.', data.expo$expo.var)] <- ell.alp.o.bet.e

    ell.alp.g.bet.o <- - t(rg) %*% (d * ro) * bet.e
    ell.alp.g.bet.o <- as.matrix(ell.alp.g.bet.o)
    hess[paste0('alp.', data.expo$geno.var), paste0('bet.', data.retro$overlap.covar)] <- ell.alp.g.bet.o

    ell.alp.g.bet.e <- ell.alp.g.bet.e - t(rg) %*% (d * prd.ro) * bet.e
    ell.alp.g.bet.e <- matrix(as.vector(ell.alp.g.bet.e), ncol = 1)
    hess[paste0('alp.', data.expo$geno.var), paste0('bet.', data.expo$expo.var)] <- ell.alp.g.bet.e

    ell.bet.o.bet.o <- - t(ro) %*% (d * ro)
    ell.bet.o.bet.o <- matrix(ell.bet.o.bet.o)
    hess[paste0('bet.', data.retro$overlap.covar), paste0('bet.', data.retro$overlap.covar)] <- ell.bet.o.bet.o

    ell.bet.o.bet.e <- - t(ro) %*% d * (alp.0 + exp(c) * bet.e) - t(ro) %*% (d * (prd.ro + prd.rg))
    ell.bet.o.bet.e <- matrix(as.vector(ell.bet.o.bet.e), ncol = 1)
    hess[paste0('bet.', data.retro$overlap.covar), paste0('bet.', data.expo$expo.var)] <- ell.bet.o.bet.e

    ell.bet.e.bet.e <- ell.bet.e.bet.e - sum(d * prd.ro) * (alp.0 + exp(c) * bet.e) - sum(d * (prd.ro^2 + 2 * prd.ro * prd.rg))
    ell.bet.e.bet.e <- matrix(ell.bet.e.bet.e)
    hess[paste0('bet.', data.expo$expo.var), paste0('bet.', data.expo$expo.var)] <- ell.bet.e.bet.e
  }

  if(length(data.retro$add.covar.retro) > 0){
    ra <- as.matrix(data.retro$rdat[, data.retro$add.covar.retro, drop = FALSE])
    ell.a.bet.a <- - t(ra) %*% d
    ell.a.bet.a <- matrix(as.vector(ell.a.bet.a), nrow = 1)
    hess['a', paste0('bet.', data.retro$add.covar.retro)] <- ell.a.bet.a

    ell.c.bet.a <- .5 * exp(c) * bet.e^2 *ell.a.bet.a
    hess['c', paste0('bet.', data.retro$add.covar.retro)] <- ell.c.bet.a

    ell.alp.0.bet.a <- bet.e * ell.a.bet.a
    hess['alp.0', paste0('bet.', data.retro$add.covar.retro)] <- ell.alp.0.bet.a

    ell.alp.g.bet.a <- - t(rg) %*% (d * ra) * bet.e
    ell.alp.g.bet.a <- matrix(ell.alp.g.bet.a)
    hess[paste0('alp.', data.expo$geno.var), paste0('bet.', data.retro$add.covar.retro)] <- ell.alp.g.bet.a

    ell.bet.a.bet.a <- - t(ra) %*% (d * ra)
    ell.bet.a.bet.a <- matrix(ell.bet.a.bet.a)
    hess[paste0('bet.', data.retro$add.covar.retro), paste0('bet.', data.retro$add.covar.retro)] <- ell.bet.a.bet.a

    ell.bet.a.bet.e <- - t(ra) %*% d * (alp.0 + exp(c) * bet.e) - t(ra) %*% (d * prd.rg)
    ell.bet.a.bet.e <- matrix(as.vector(ell.bet.a.bet.e), ncol = 1)
    hess[paste0('bet.', data.retro$add.covar.retro), paste0('bet.', data.expo$expo.var)] <- ell.bet.a.bet.e
  }

  if(length(data.retro$overlap.covar) > 0 && length(data.retro$add.covar.retro) > 0){
    ell.alp.o.bet.a <- - t(ro) %*% (d * ra) * bet.e
    ell.alp.o.bet.a <- matrix(ell.alp.o.bet.a)
    hess[paste0('alp.', data.expo$overlap.covar), paste0('bet.', data.retro$add.covar.retro)] <- ell.alp.o.bet.a

    ell.bet.a.bet.o <- - t(ra) %*% (d * ro)
    ell.bet.a.bet.o <- matrix(ell.bet.a.bet.o)
    hess[paste0('bet.', data.retro$add.covar.retro), paste0('bet.', data.retro$overlap.covar)] <- ell.bet.a.bet.o

    ell.bet.a.bet.e <- ell.bet.a.bet.e - t(ra) %*% (d * prd.ro)
    ell.bet.a.bet.e <- matrix(as.vector(ell.bet.a.bet.e), ncol = 1)
    hess[paste0('bet.', data.retro$add.covar.retro), paste0('bet.', data.expo$expo.var)] <- ell.bet.a.bet.e
  }

  ell.c.bet.e <- exp(c) * bet.e * ell.a + .5 * exp(c) * bet.e^2 * ell.a.bet.e
  hess['c', paste0('bet.', data.expo$expo.var)] <- ell.c.bet.e

  ell.alp.0.bet.e <- ell.a + bet.e * ell.a.bet.e
  ell.alp.0.bet.e <- matrix(ell.alp.0.bet.e)
  hess['alp.0', paste0('bet.', data.expo$expo.var)] <- ell.alp.0.bet.e

  ell.bet.e.bet.e <- ell.bet.e.bet.e + ell.a.bet.e * (alp.0 + exp(c) * bet.e)
  ell.bet.e.bet.e <- matrix(ell.bet.e.bet.e)
  hess[paste0('bet.', data.expo$expo.var), paste0('bet.', data.expo$expo.var)] <- ell.bet.e.bet.e

  for(i in 1:nrow(hess)){
    for(j in 1:ncol(hess)){
      if(is.na(hess[i, j])){
        hess[i, j] <- hess[j, i]
      }else{
        break
      }
    }
  }

  hess <- -hess

  list(gr = gr, hess = hess)

}

