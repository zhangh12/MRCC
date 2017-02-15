
check.LMT.ci <- function(rdata, edata, c.lmt = 1.0, level, b.ci){

  if(is.null(b.ci)){
    return(NA)
  }

  crit <- qchisq(level, df = 1)
  par.pos <- align.parameter(find.null.tsls(rdata, edata))$par.pos
  name.bet.z <- paste0('bet.', edata$vz)

  if(abs(b.ci) < 1e-3){
    par.mle <- align.parameter(find.null.tsls(rdata, edata))$par
  }else{
    par.mle <-  try(find.LRT.mle(rdata, edata, b.ci), silent = TRUE)
    if('try-error' %in% class(par.mle)){
      return(NA)
    }
    par.mle <- c(par.mle, b.ci)
    names(par.mle)[par.pos[name.bet.z, 'start']] <- name.bet.z
  }

  J0 <- -hessian0(par.mle, rdata, edata, par.pos)
  I0 <- fisher.info0(par.mle, rdata, edata, par.pos)
  sc <- score(par.mle, rdata, edata, par.pos)[name.bet.z]

  id.bet <- which(colnames(J0) == name.bet.z)
  v <- cbind(-J0[id.bet,-id.bet] %*% solve(J0[-id.bet, -id.bet]), 1)
  colnames(v)[ncol(v)] <- name.bet.z
  v <- v[, colnames(I0), drop = FALSE]

  lm.stat <- sc^2/(v %*% I0 %*% t(v))

  ifelse(lm.stat <= crit, TRUE, FALSE)

}

