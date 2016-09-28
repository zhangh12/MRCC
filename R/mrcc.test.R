
mrcc.test <- function(rdata, edata, tpar = NULL){

  nm <- naive.method(rdata, edata)
  naive.est <- nm$naive.est
  summary.h <- nm$summary.h
  edata <- nm$edata # vz might have been changed as residuals

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
  #   print(c(range(h1-h2)))
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
  #   print(range(tmp/gr)-1)
  # }

  fit <- newton.raphson(par, rdata, edata, par.pos)

  foo <- function(par, rdata){

    alp.g <- par[paste0('alp.', rdata$vg)]
    rg <- as.matrix(rdata$data[, rdata$vg, drop = FALSE])
    u <- rg %*% alp.g
    d <- rdata$data[, rdata$vd, drop = FALSE]
    dat <- data.frame(d, u)
    fit <- glm(paste(rdata$vd,'~ u'), data = dat, family = 'binomial')
    print(coef(fit)['u'] - par[paste0('bet.', edata$vz)])

  }

  #foo(par, rdata)

  if(!fit$convergence){
    return(list(res = NULL, gr = NULL, fn = NULL, convergence = fit$convergence))
  }

  hess <- hessian(fit$par, rdata, edata, par.pos)
  hess0 <- hessian0(fit$par, rdata, edata, par.pos)
  info <- fisher.info(fit$par, rdata, edata, par.pos)
  info0 <- fisher.info0(fit$par, rdata, edata, par.pos)

  inv.hess <- solve(hess)
  cov1 <- inv.hess %*% info %*% inv.hess
  se1 <- sqrt(diag(cov1))
  names(se1) <- names(fit$par)

  inv.hess0 <- solve(hess0)
  cov0 <- inv.hess0 %*% info0 %*% inv.hess0
  se <- sqrt(diag(cov0))
  names(se) <- names(fit$par)

  if(!is.null(tpar)){
    thess <- hessian0(tpar, rdata, edata, par.pos)
    tinfo <- fisher.info0(tpar, rdata, edata, par.pos)
    inv.thess <- solve(thess)
    tcov <- inv.thess %*% tinfo %*% inv.thess
    tse <- sqrt(diag(tcov))
    names(tse) <- names(fit$par)
  }

  if(is.null(tpar)){
    summary.n <- data.frame(Estimate = fit$par, SE = se, SE1 = NA, stringsAsFactors = FALSE)
  }else{
    summary.n <- data.frame(Estimate = fit$par, SE = se, SE1 = tse, stringsAsFactors = FALSE)
  }

  rownames(summary.n) <- names(fit$par)
  summary.n$z <- summary.n$Estimate / summary.n$SE
  summary.n$"Pr(>|z|)" <- pchisq(summary.n$z^2, df = 1, lower.tail = FALSE)

  res <- rbind(summary.n, summary.h)

  var.e <- paste0('alp.', c('0', edata$vh, edata$vx, edata$vg))
  var.r <- paste0('bet.', c(rdata$vx, rdata$vy, edata$vz))
  var.0 <- setdiff(names(par), c(var.e, var.r))
  #var.0 <- NULL

  res <- res[c(var.0, var.e, var.r), ]


  #####
  p.lrt <- LRT(rdata, edata, fit$fn)

  list(par = fit$par, MRCC = res, p.lrt = p.lrt, gr = fit$gr, fn = fit$fn, convergence = fit$convergence)

}

