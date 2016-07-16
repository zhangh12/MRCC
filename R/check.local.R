
check.local <- function(par, rdata, edata, par.pos){

  s <- seq(-2,2,length.out = 100)
  par(mfrow=c(4,4))
  for(kk in 1:length(par)){
    f <- NULL
    for(s0 in s){
      par1 <- par
      par1[kk] <- s0
      f <- c(f, neg.log.lik(par1, rdata, edata, par.pos))
    }
    plot(s,f,type='p',pch=20,cex=.1,xlab = names(par)[kk])
    abline(v=par[kk],col='red')
  }

}