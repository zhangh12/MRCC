
check.mle <- function(par, rdata, edata, par.pos, tol = 1e-8){

  if(any(is.na(par))){
    return(FALSE)
  }

  t1 <- try(s <- score(par, rdata, edata, par.pos))
  if('try-error' %in% class(t1)){
    return(FALSE)
  }

  if(max(abs(s)) > tol){
    return(FALSE)
  }

  t2 <- try(inv.h <- solve(hessian(par, rdata, edata, par.pos)))
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
