
LM.test <- function(rdata, edata, par, se, level, b.ci, fig){

  null.mle <- find.null.tsls(rdata, edata)
  ap <- align.parameter(null.mle)
  par.mle <- ap$par
  par.pos <- ap$par.pos

  name.bet.z <- paste0('bet.', edata$vz)

  J0 <- -hessian0(par.mle, rdata, edata, par.pos)
  #J <- -hessian(par.mle, rdata, edata, par.pos)
  sc <- score(par.mle, rdata, edata, par.pos)[name.bet.z]

  if(is.null(rdata$vx) & is.null(rdata$vy)){
    id.bet <- which(colnames(J0) == name.bet.z)
    lm.stat <- sc^2 * solve(J0)[id.bet, id.bet]
    #lm.stat <- sc^2 * solve(J)[id.bet, id.bet]
    p.lmt <- pchisq(lm.stat, df = 1, lower.tail = FALSE)
  }else{

  }

  ci <- LMT.conf.int(rdata, edata, par, se, c.lmt = 1.0, level, fig)

  in.ci <- check.LMT.ci(rdata, edata, c.lmt = 1.0, level, b.ci)

  list(p.lmt = p.lmt, ci = ci, in.ci = in.ci, stat = lm.stat)

}
