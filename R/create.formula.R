

create.formula <- function(rdata, edata){


  rform <- paste(rdata$vd, '~ pred.expo')
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
  eform <- as.formula(eform)

  if(length(edata$vh) > 0){
    nform <- paste(edata$vz, '~', paste(edata$vh, collapse = '+'))
    nform <- as.formula(nform)
  }else{
    nform <- NULL
  }

  list(rform = rform, eform = eform, nform = nform)

}
