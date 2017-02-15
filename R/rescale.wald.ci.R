
rescale.wald.ci <- function(res.wald, res.lm, edata, level, b.ci){

  ci1 <- res.wald$ci
  in.ci1 <- res.wald$in.ci

  wald.stat <- res.wald$stat
  lm.stat <- res.lm$stat

  c.wald <- wald.stat/lm.stat

  name.bet.z <- paste0('bet.', edata$vz)
  ci <- res.wald$par[name.bet.z] + c(-1, 1) * res.wald$se[name.bet.z] * sqrt(c.wald * qchisq(level, df = 1))
  names(ci) <- c('LCL', 'RCL')

  in.ci <- ifelse(is.null(b.ci), NA, ci['LCL'] <= b.ci && ci['RCL'] >= b.ci)

  res.wald$c.wald <- c.wald
  res.wald$ci <- ci
  res.wald$in.ci <- in.ci
  res.wald$ci1 <- ci1
  res.wald$in.ci1 <- in.ci1

  res.wald

}
