

create.null.formula <- function(rdata, edata){

  if(length(rdata$vx) == 0 & length(rdata$vy) == 0){
    rform <- paste(rdata$vd, '~ 1')
  }else{
    if(length(rdata$vx) != 0 & length(rdata$vy) == 0){
      rform <- paste(rdata$vd, '~ ', paste(rdata$vx, collapse = '+'))
    }else{
      if(length(rdata$vx) == 0 & length(rdata$vy) != 0){
        rform <- paste(rdata$vd, '~ ', paste(rdata$vy, collapse = '+'))
      }else{
        rform <- paste(rdata$vd, '~ ', paste(c(rdata$vx, rdata$vy), collapse = '+'))
      }
    }
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
