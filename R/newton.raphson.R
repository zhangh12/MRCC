
newton.raphson <- function(par, rdata, edata, par.pos, v = 0, tol = 1e-8, max.it = 200){

  fn <- logL(par, rdata, edata, par.pos)
  it <- 0
  while(TRUE){

    it <- it + 1
    if(it > max.it){
      stop('Newton-Raphson algorithm fails')
    }

    h <- hessian(par, rdata, edata, par.pos)
    s <- score(par, rdata, edata, par.pos, v)
    tau <- as.vector(solve(h) %*% s)
    names(tau) <- names(s)
    ncut <- 20
    cut <- seq(0, 1, length.out = ncut+1)[-1]

    e <- NULL
    f <- NULL
    for(i in 1:ncut){
      par2 <- par - cut[i] * tau
      s2 <- score(par2, rdata, edata, par.pos, v)
      e <- c(e, max(abs(s2)))
      f <- c(f, logL(par2, rdata, edata, par.pos))
    }

    if(is.null(e)){
      stop('Newton-Raphson algorithm fails')
    }

    if(length(fn) > 20 & max(fn) - max(f) > 1e-2){
      stop('Newton-Raphson algorithm fails')
    }

    if(all(v == 0)){
      par <- par - cut[which.max(f)] * tau
    }else{
      par <- par - cut[which.min(e)] * tau
    }
    fn <- c(fn, logL(par, rdata, edata, par.pos))
    #plot(fn, type = 'b', pch = 20, cex = .2)
    s <- score(par, rdata, edata, par.pos, v)
    eps <- max(abs(s))
    #print(eps)
    if(eps < tol){
      #message('converged')
      break
    }
  }

  par

}