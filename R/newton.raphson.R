
newton.raphson <- function(par, rdata, edata, par.pos, tol = 1e-8, max.it = 100, max.seed = 10){

  set.seed(0)
  par0 <- par
  eps <- 1
  it <- 0
  nseed <- 0
  convergence <- 0
  f0 <- neg.log.lik(par0, rdata, edata, par.pos)
  f <- f0

  par.c <- NULL # candidate
  fn.c <- NULL
  while(nseed <= max.seed && !convergence){
    while(TRUE){
      hess <- hessian(par, rdata, edata, par.pos)

      if('try-error' %in% class(try(inv.hess <- solve(hess), silent = TRUE))){
        message('cannot apply newton\'s method when inversing hessian')
        par <- reset.ini(par0, nseed)
        #readline('@')
        eps <- 1
        it <- 0
        nseed <- nseed + 1
        convergence <- 0
        f <- f0
        break
      }

      gr <- score(par, rdata, edata, par.pos)
      delta <- as.vector(inv.hess %*% gr)

      par1 <- par
      ls <- line.search(par, rdata, edata, par.pos, delta)
      status.ls <- ls$status
      par <- ls$par
      cd <- coordinate.descent(par, rdata, edata, par.pos)
      status.cd <- cd$status
      par <- cd$par

      eps.par <- max(abs(par1 - par))
      eps.gr <- max(abs(score(par, rdata, edata, par.pos)))
      eps <- max(eps.par, eps.gr)

      #print(c(eps.par, eps.gr))
      if(eps <= 1e-4){
        if(valid.solution(par, rdata, edata, par.pos)){
          #message('keep a candidate')
          if(is.null(par.c)){
            par.c <- par
            fn.c <- neg.log.lik(par, rdata, edata, par.pos)
          }else{
            tmp <- neg.log.lik(par, rdata, edata, par.pos)
            if(tmp < fn.c){
              par.c <- par
              fn.c <- tmp
            }
          }
        }
      }

      # print(par1-par)
      # print(score(par, rdata, edata, par.pos))
      # print(neg.log.lik(par, rdata, edata, par.pos))
      #check.local(par, rdata, edata, par.pos)
      ##readline('@')

      if(status.ls || status.cd){
        f <- c(f, neg.log.lik(par, rdata, edata, par.pos))
      }else{
        #message('both line search and coordinate-wise search fail')

        if(!is.null(par.c)){
          convergence <- 2
          break
        }

        par <- reset.ini(par0, nseed)
        #readline('@')
        eps <- 1
        it <- 0
        nseed <- nseed + 1
        convergence <- 0
        f <- f0
        break
      }

      #plot(f, type = 'b', pch = 4, ylim = range(f), xlab = 'iter', ylab = '-log L', col = 'red', cex = .3)

      #readline('>>')

      if(eps <= tol){ # is it a real solution?
        conv <- valid.solution(par, rdata, edata, par.pos)
        if(!conv){ # no
          message('not a valid solution even if gradient is zero')
          par <- reset.ini(par0, nseed)
          #readline('@')
          eps <- 1
          it <- 0
          nseed <- nseed + 1
          convergence <- 0
          f <- f0
          break
        }else{ # yes
          convergence <- 1
          break
        }
      }

      it <- it + 1

      if(it > max.it){
        message('maximum interation step achieved')
        par <- reset.ini(par0, nseed)
        #readline('@')
        eps <- 1
        it <- 0
        nseed <- nseed + 1
        convergence <- 0
        f <- f0
        break
      }
    }
  }

  if(!convergence){
    if(is.null(par.c)){
      par[] <- NA
    }else{
      par <- par.c
    }
  }

  #check.local(par, rdata, edata, par.pos)
  #plot(f, type = 'b', pch = 4, ylim = range(f), xlab = 'iter', ylab = '-log L', col = 'red', cex = .3)
  ##readline('@')

  fn <- neg.log.lik(par, rdata, edata, par.pos)
  gr <- score(par, rdata, edata, par.pos)
  hess <- hessian(par, rdata, edata, par.pos)
  # dhess <- optim(par, neg.log.lik, gr = score,
  #                rdata = rdata, edata = edata, par.pos = par.pos,
  #                method = 'BFGS', control = list(trace = 0, maxit = 0),
  #                hessian = TRUE)$hessian

  list(par = par, gr = gr, hess = hess, fn = fn, convergence = convergence)

}