
# find MLE by rotation, bet -> alp -> bet -> ...
find.mle.2S <- function(par, rdata, edata, par.pos){

  bet <- par[paste0('bet.', edata$vz)]
  rg <- as.matrix(rdata$rg)
  n <- nrow(rg)
  n1 <- sum(rdata$rd)
  n0 <- n - n1

  f <- logL(par, rdata, edata, par.pos)

  while(TRUE){

    par <- find.LRT.mle(rdata, edata, bet) # update par except for bet
    alp <- par[paste0('alp.', rdata$vg)]
    rdata$data$o <- log(n1/n0) + par['a']
    rdata$data$pred.expo <- rg %*% alp
    sform <- paste(rdata$vd, '~ offset(o) + pred.expo - 1')
    sfit <- glm(sform, data = rdata$data, family = 'binomial')
    bet <- coef(sfit)['pred.expo']
    s <- score.LRT(par, rdata, edata, par.pos, bet)
    f <- c(f, logL.LRT(par, rdata, edata, par.pos, bet))
    plot(f, type = 'b')
    print(s)
    if(all(abs(s) < 1e-7)){
      break
    }

  }

  par

}

