
logL.ind <- function(par, rdata, edata, par.pos){

  if(any(is.na(par))){
    return(NULL)
  }

  ne <- nrow(edata$data)
  nr <- nrow(rdata$data)

  a <- par['a']

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

  lin <- a

  rg <- as.matrix(rdata$data[, rdata$vg, drop = FALSE])
  eg <- as.matrix(edata$data[, edata$vg, drop = FALSE])

  lin <- lin + rg %*% alp.g * bet.z
  if(length(rdata$vx) > 0){
    rx <- as.matrix(rdata$data[, rdata$vx, drop = FALSE])
    lin <- lin + rx %*% bet.x
  }

  if(length(rdata$vy) > 0){
    ry <- as.matrix(rdata$data[, rdata$vy, drop = FALSE])
    lin <- lin + ry %*% bet.y
  }

  lin <- as.vector(lin)
  delta <- exp(lin)

  rd <- as.vector(rdata$data[, rdata$vd])
  n1 <- sum(rd)
  n0 <- nr - n1

  r <- edata$data[, edata$vz] - alp.0 - eg %*% alp.g

  if(length(edata$vx) > 0){
    ex <- as.matrix(edata$data[, edata$vx, drop = FALSE])
    r <- r - ex %*% alp.x
  }
  r <- as.vector(r)

  ld <- rd * lin - log(n0) - log(1 + n1/n0 * delta)
  ld1 <- ld[rd == 1]
  ld0 <- ld[rd == 0]
  lz <- -log(2*pi)/2 - log(c)/2 - 1/c/2 * r^2
  ll <- data.frame(val=c(ld1,ld0,lz),
                   grp = rep(c('cc1', 'cc0','exp'),times=c(length(ld1), length(ld0),length(lz))))
  #boxplot(val~grp,data=ll)

  l <- sum(rd * lin) - nr * log(n0) - sum(log(1 + n1/n0 * delta)) - ne * log(2*pi)/2 - ne * log(c)/2 - 1/c/2 * sum(r^2)
  names(l) <- NULL
  l

}


