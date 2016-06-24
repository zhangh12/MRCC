
mrcc.test <- function(data.retro, data.expo, naive.est){

  ap <- align.parameter(naive.est)
  par <- ap$par
  par.pos <- ap$par.pos

  # for(i in 1:1000){
  #   par1 <- par + runif(length(par), -.2, .2)
  #   par2 <- par1
  #   par2['bet.z'] <- par2['bet.z'] + 1e-8
  #   gr1 <- score.hess(par1, data.retro, data.expo, par.pos)$gr
  #   gr2 <- score.hess(par2, data.retro, data.expo, par.pos)$gr
  #   h1 <- (gr2['bet.z']-gr1['bet.z'])/1e-8
  #   h2 <- score.hess(par1, data.retro, data.expo, par.pos)$hess['bet.z', 'bet.z']
  #   print(h1/h2)
  # }

  # for(i in 1:1000){
  #   par1 <- par + runif(length(par), -.5, .5)
  #   gr <- score.hess(par1, data.retro, data.expo, par.pos)$gr
  #   tmp <- NULL
  #   for(j in 1:length(par)){
  #     par2 <- par1
  #     par2[j] <- par2[j] + 1e-8
  #     fn1 <- neg.log.lik(par1, data.retro, data.expo, par.pos)
  #     fn2 <- neg.log.lik(par2, data.retro, data.expo, par.pos)
  #     tmp <- c(tmp, (fn2-fn1)/1e-8)
  #   }
  #   print(tmp/gr)
  #   print('')
  # }

  if(0){
    fit0 <- optim(par, neg.log.lik, gr = NULL, data.retro = data.retro, data.expo = data.expo, par.pos = par.pos, method = 'BFGS', control = list(trace = 0, maxit = 1e4), hessian = TRUE)
    se0 <- sqrt(diag(solve(fit0$hessian)))
    names(se0) <- names(fit0$par)
  }

  fit <- newton.raphson(par, data.retro, data.expo, par.pos)

  info <- solve(fit$hessian)
  se <- sqrt(diag(info))
  names(se) <- names(fit$par)

  summary.n <- data.frame(Estimate = fit$par, SE = se, stringsAsFactors = FALSE)
  rownames(summary.n) <- names(par)
  summary.n$z <- summary.n$Estimate / summary.n$SE
  summary.n$"Pr(>|z|)" <- pchisq(summary.n$z^2, df = 1, lower.tail = FALSE)

  list(summary.n = summary.n, gr = fit$gr, fn = fit$fn, convergence = fit$convergence)

}

