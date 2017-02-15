
# check whether a solution of MLE is valid
# to be a valid MLE, it should be
# 1. no NA
# 2. has scores close to zero
# 3. invertible hessian matrix
# 4. positive estimates of variance for all parameters
check.mle <- function(par, rdata, edata, par.pos, tol = 1e-8){

  if(is.null(par) || any(is.na(par))){
    return(FALSE)
  }

  if(par['c'] <= 0){
    return(FALSE)
  }

  t1 <- try(s <- score(par, rdata, edata, par.pos))
  if('try-error' %in% class(t1)){
    return(FALSE)
  }

  # if(max(abs(s)) > tol){
  #   return(FALSE)
  # }

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
