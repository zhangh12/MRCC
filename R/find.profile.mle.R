

find.profile.mle <- function(rdata, edata){

  bet <- seq(-2, 2, by = .05)
  b0 <- NULL
  fn <- NULL
  for(b in bet){
    ini <- find.LRT.tsls(rdata, edata, b)

    ap <- align.parameter(ini, null = TRUE)
    par.ini <- ap$par
    par.pos <- ap$par.pos

    t1 <- try(sol <- optimx(par.ini, nlogL.LRT, gr = nscore.LRT, method = 'ucminf',
                            rdata = rdata, edata = edata, par.pos = par.pos, bet = b), silent = TRUE)

    # t1 <- try(sol <- optimx(par.ini, nlogL.LRT, gr = nscore.LRT, hess = nhessian.LRT, method = 'ucminf',
    #                         rdata = rdata, edata = edata, par.pos = par.pos, bet = b), silent = TRUE)

    if('try-error' %in% class(t1)){
      next
    }

    if(!sol['convcode']){
      par <- as.vector(t(sol[1, names(par.ini)]))
      names(par) <- names(par.ini)

      if(check.LRT.mle(par, rdata, edata, par.pos, b)){
        b0 <- c(b0, b)
        fn <- c(fn, logL.LRT(par, rdata, edata, par.pos, b))
      }
    }

    plot(b0, fn, type = 'b', pch = 20, cex = .3, main = 'Profile Likelihood of beta')
  }

  b0 <- b0[which.max(fn)]
  abline(v = b0)
  par <- find.LRT.mle(rdata, edata, b0)
  par[paste0('bet.', edata$vz)] <- b0
  par <- find.mle(par, rdata, edata, align.parameter(find.tsls(rdata, edata))$par.pos)

  par

}
