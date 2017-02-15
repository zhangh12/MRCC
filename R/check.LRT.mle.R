
check.LRT.mle <- function(par, rdata, edata, par.pos, bet, tol = 1e-8){

  if(any(is.na(par))){
    return(FALSE)
  }

  if(par['c'] <= 0){
    return(FALSE)
  }

  t1 <- try(s <- score.LRT(par, rdata, edata, par.pos, bet))
  if('try-error' %in% class(t1)){
    return(FALSE)
  }

  # if(max(abs(s)) > tol){
  #   return(FALSE)
  # }

  t2 <- try(inv.h <- solve(hessian.LRT(par, rdata, edata, par.pos, bet)))
  if('try-error' %in% class(t2)){
    return(FALSE)
  }

  if(any(is.na(inv.h))){
    return(FALSE)
  }

  se <- -diag(inv.h)
  if(any(se <= 0)){
    return(FALSE)
  }

  return(TRUE)

}
