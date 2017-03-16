

create.LRT.formula <- function(rdata, edata, bet){

  rform <- paste(rdata$vd, '~ offset(I(', bet, '* pred.expo))')
  if(length(rdata$vx) > 0){
    rform <- paste(rform, '+', paste(rdata$vx, collapse = '+'))
  }
  if(length(rdata$vy) > 0){
    rform <- paste(rform, '+', paste(rdata$vy, collapse = '+'))
  }
  rform <- as.formula(rform)

  eform <- paste(edata$vz, '~', edata$vd, '+', paste(edata$vg, collapse = '+'))
  if(length(edata$vx) > 0){
    eform <- paste(eform, '+', paste(edata$vx, collapse = '+'))
  }
  if(length(edata$vh) > 0){
    eform <- paste(eform, '+', paste(edata$vh, collapse = '+'))
  }
  eform <- as.formula(eform)

  list(rform = rform, eform = eform)

}
