
newton.raphson.LRT <- function(par, rdata, edata, par.pos, v = 0, tol = 1e-8){

  id <- which(names(par) == paste0('bet.', edata$vz))

  while(TRUE){

    h <- hessian.LRT(par, rdata, edata, par.pos)
    s <- score.LRT(par, rdata, edata, par.pos, v)
    tau <- as.vector(solve(h) %*% s)
    names(tau) <- names(s)
    ncut <- 20
    cut <- seq(0, 1, length.out = ncut+1)[-1]

    e <- NULL
    for(i in 1:ncut){
      par2 <- par
      par2[-id] <- par[-id] - cut[i] * tau
      s2 <- score.LRT(par2, rdata, edata, par.pos, v)
      e <- c(e, max(abs(s2)))
    }

    if(is.null(e)){
      stop('Newton-Raphson algorithm fails')
    }else{
      par[-id] <- par[-id] - cut[which.min(e)] * tau
      s <- score.LRT(par, rdata, edata, par.pos, v)
      eps <- max(abs(s))
      #print(eps)
      if(eps < tol){
        #message('converged')
        break
      }
    }
  }

  par

}
