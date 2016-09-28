
TSLS.test <- function(rdata, edata, tpar = NULL){

  tmp <- TSLS(rdata, edata)
  naive.est <- tmp$naive.est
  summary.h <- tmp$summary.h
  edata <- tmp$edata

  ap <- align.parameter(naive.est)
  par <- ap$par
  par.pos <- ap$par.pos

  hess <- TSLS.hessian(par, rdata, edata, par.pos)
  info <- TSLS.fisher.info(par, rdata, edata, par.pos)

  cov <- solve(hess) %*% info %*% t(solve(hess))
  se <- sqrt(diag(cov))
  names(se) <- names(par)

  if(!is.null(tpar)){
    thess <- TSLS.hessian(tpar, rdata, edata, par.pos)
    tinfo <- TSLS.fisher.info(tpar, rdata, edata, par.pos)
    inv.thess <- solve(thess)
    tcov <- inv.thess %*% tinfo %*% inv.thess
    tse <- sqrt(diag(tcov))
    names(tse) <- names(par)
  }

  if(is.null(tpar)){
    summary.n <- data.frame(Estimate = par, SE = se, SE1 = NA, stringsAsFactors = FALSE)
  }else{
    summary.n <- data.frame(Estimate = par, SE = se, SE1 = tse, stringsAsFactors = FALSE)
  }

  rownames(summary.n) <- names(par)
  summary.n$z <- summary.n$Estimate / summary.n$SE
  summary.n$"Pr(>|z|)" <- pchisq(summary.n$z^2, df = 1, lower.tail = FALSE)

  res <- rbind(summary.n, summary.h)

  var.e <- paste0('alp.', c('0', edata$vh, edata$vx, edata$vg))
  var.r <- paste0('bet.', c(rdata$vx, rdata$vy, edata$vz))
  var.0 <- setdiff(names(par), c(var.e, var.r))
  #var.0 <- NULL

  res <- res[c(var.0, var.e, var.r), ]

  res

}
