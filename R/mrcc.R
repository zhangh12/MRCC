

mrcc <- function(rformula, rdata, eformula, edata){

  data.retro <- data.parse.retro(rformula, rdata)
  data.expo <- data.parse.expo(eformula, edata)

  exd <- extend.data(data.retro, data.expo)

  data.retro <- exd$data.retro
  data.expo <- exd$data.expo

  nm <- naive.method(data.retro, data.expo)
  naive.est <- nm$naive.est
  summary.a <- nm$summary.a

  mt <- mrcc.test(data.retro, data.expo, naive.est)
  summary.n <- mt$summary.n
  convergence <- mt$convergence

  res <- rbind(summary.n, summary.a)

  var.e <- c('alp.0', paste0('alp.', data.expo$covar.var), paste0('alp.', data.expo$geno.var))
  var.r <- c(paste0('bet.', data.retro$covar.var), paste0('bet.', data.expo$expo.var))

  var.e <- intersect(var.e, rownames(res))
  var.r <- intersect(var.r, rownames(res))
  var.0 <- setdiff(rownames(res), c(var.e, var.r))

  res <- res[c(var.0, var.e, var.r), ]

  list(res = res, gr = mt$gr, fn = mt$fn, convergence = convergence)

}

