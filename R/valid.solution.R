
valid.solution <- function(par, rdata, edata, par.pos){

  hess <- hessian(par, rdata, edata, par.pos)
  info <- fisher.info(par, rdata, edata, par.pos)
  #info0 <- fisher.info0(par, rdata, edata, par.pos)
  #ohess <- optim(par, neg.log.lik, gr = score,rdata = rdata, edata = edata, par.pos = par.pos,method = 'BFGS', control = list(trace = 0, maxit = 0),hessian = TRUE)$hessian
  tt <- try(inv.hess <- solve(hess), silent = TRUE)
  ss <- try(inv.info <- solve(info), silent = TRUE)

  ifelse(any(is.na(hess)) || ('try-error' %in% class(tt)) || ('try-error' %in% class(ss)) || any(is.na(inv.hess)) || any(diag(inv.hess) < 0) || any(is.na(inv.info)) || any(diag(inv.info) < 0), FALSE, TRUE)

}
