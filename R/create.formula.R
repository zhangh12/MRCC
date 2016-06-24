

create.formula <- function(data.retro, data.expo){


  rform <- paste(data.retro$retro.var, '~ pred.expo')
  if(length(data.retro$covar.var) > 0){
    rform <- paste(rform, '+', paste(data.retro$covar.var, collapse = '+'))
  }
  rform <- as.formula(rform)

  eform <- paste(data.expo$expo.var, '~', data.expo$retro.var, '+', paste(data.expo$geno.var, collapse = '+'))
  if(length(data.expo$overlap.covar) > 0){
    eform <- paste(eform, '+', paste(data.expo$overlap.covar, collapse = '+'))
  }
  eform <- as.formula(eform)

  if(length(data.expo$add.covar.expo) > 0){
    nform <- paste(data.expo$expo.var, '~', paste(data.expo$add.covar.expo, collapse = '+'), '-1')
    nform <- as.formula(nform)
  }else{
    nform <- NULL
  }

  list(rform = rform, eform = eform, nform = nform)

}
