
neg.log.lik <- function(par, data.retro, data.expo, par.pos){

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
  if(max(lin) > 10){
    offset <- max(lin) - 10
  }else{
    offset <- 0
  }

  lin <- lin - offset

  delta <- exp(lin)

  n1 <- sum(data.retro$rdat[, data.retro$retro.var])
  n0 <- nr - n1
  p <- (exp(-offset) + n1/n0 * delta)

  res <- data.expo$edat[, data.expo$expo.var] - alp.0 - as.matrix(data.expo$edat[, data.expo$geno.var]) %*% alp.g
  if(!is.na(b)){
    res <- res - (bet.e * exp(c) + b) * data.expo$edat[, data.expo$retro.var]
  }

  if(length(data.expo$overlap.covar) > 0){
    res <- res - as.matrix(data.expo$edat[, data.expo$overlap.covar, drop = FALSE]) %*% alp.o
  }

  -(-nr*(offset+log(n0)) - sum(log(p)) + sum(data.retro$rdat[, data.retro$retro.var] * lin) - .5 * log(2 * pi) * ne - .5 * c * ne -
      .5 / exp(c) * sum(res^2))

}


