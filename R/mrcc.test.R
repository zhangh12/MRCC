
mrcc.test <- function(rdata, edata){

  nm <- naive.method(rdata, edata)
  naive.est <- nm$naive.est
  summary.h <- nm$summary.h
  edata <- nm$edata # vz might have been changed

  ap <- align.parameter(naive.est)
  par <- ap$par
  par.pos <- ap$par.pos

  # for(i in 1:1000){
  #   par1 <- par + runif(length(par), -.2, .2)
  #   h1 <- optim(par1, neg.log.lik, gr = score,
  #               rdata = rdata, edata = edata, par.pos = par.pos,
  #               method = 'BFGS', control = list(trace = 0, maxit = 0),
  #               hessian = TRUE)$hessian
  #   h2 <- hessian(par1, rdata, edata, par.pos)
  #   print(c(i,range(h1/h2)))
  # }

  # for(i in 1:10000){
  #   par1 <- par + runif(length(par), -1, 1)
  #   gr <- score(par1, rdata, edata, par.pos)
  #   tmp <- NULL
  #   for(j in 1:length(par)){
  #     par2 <- par1
  #     par2[j] <- par2[j] + 1e-8
  #     fn1 <- neg.log.lik(par1, rdata, edata, par.pos)
  #     fn2 <- neg.log.lik(par2, rdata, edata, par.pos)
  #     tmp <- c(tmp, (fn2-fn1)/1e-8)
  #   }
  #   print(range(tmp/gr))
  # }

  fit <- newton.raphson(par, rdata, edata, par.pos)

  if(!fit$convergence){
    return(list(res = NULL, gr = NULL, fn = NULL, convergence = fit$convergence))
  }

  hess <- fit$hess
  info <- fisher.info(fit$par, rdata, edata, par.pos)

  inv.hess <- solve(hess)
  cov <- inv.hess %*% info %*% inv.hess
  se <- sqrt(diag(cov))
  names(se) <- names(fit$par)

  se0 <- sqrt(diag(inv.hess))
  names(se0) <- names(fit$par)

  summary.n <- data.frame(Estimate = fit$par, SE = se, SE0 = se0, stringsAsFactors = FALSE)
  rownames(summary.n) <- names(par)
  summary.n$z <- summary.n$Estimate / summary.n$SE
  summary.n$"Pr(>|z|)" <- pchisq(summary.n$z^2, df = 1, lower.tail = FALSE)

  res <- rbind(summary.n, summary.h)

  var.e <- paste0('alp.', c('0', edata$vh, edata$vx, edata$vg))
  var.r <- paste0('bet.', c(rdata$vx, rdata$vy, edata$vz))
  var.0 <- setdiff(names(par), c(var.e, var.r))

  res <- res[c(var.0, var.e, var.r), ]

  list(res = res, gr = fit$gr, fn = fit$fn, convergence = fit$convergence)

}

