
newton.raphson <- function(par, data.retro, data.expo, par.pos, tol = 1e-8, max.it = 1e5, max.seed = 100){

  set.seed(0)
  par0 <- par
  eps <- 1
  it <- 0
  nseed <- 0
  convergence <- 0
  while(nseed <= max.seed && !convergence){
    while(TRUE){
      sh <- score.hess(par, data.retro, data.expo, par.pos)

      if('try-error' %in% class(try(inv <- solve(sh$hess), silent = TRUE))){
        par <- reset.ini(par0, nseed)
        eps <- 1
        it <- 0
        nseed <- nseed + 1
        convergence <- 0
        break
      }

      delta <- as.vector(inv %*% sh$gr)
      eps.par <- max(abs(delta))
      eps.gr <- max(abs(sh$gr))
      eps <- max(eps.par, eps.gr)

      par[] <- par - .1 * delta
      it <- it + 1

      if(it > max.it){
        par <- reset.ini(par0, nseed)
        eps <- 1
        it <- 0
        nseed <- nseed + 1
        convergence <- 0
        break
      }

      if(eps <= tol){ # is it a real solution?
        sh <- score.hess(par, data.retro, data.expo, par.pos)
        tt <- try(inv <- solve(sh$hess), silent = TRUE)
        if(any(is.na(sh$hess)) || ('try-error' %in% class(tt)) || any(is.na(inv)) || any(diag(inv) < 0)){ # no
          par <- reset.ini(par0, nseed)
          eps <- 1
          it <- 0
          nseed <- nseed + 1
          convergence <- 0
          break
        }else{ # yes
          fn <- neg.log.lik(par, data.retro, data.expo, par.pos)
          convergence <- 1
          break
        }
      }
    }
  }

  if(!convergence){
    par[] <- NA
    gr <- NULL
    hessian <- NULL
    fn <- NULL
  }else{
    gr <- sh$gr
    hess <- sh$hess
  }

  list(par = par, gr = gr, hessian = hess, fn = fn, convergence = convergence)

}