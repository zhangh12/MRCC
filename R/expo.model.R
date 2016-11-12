
expo.model <- function(edata){

  eform <- paste(edata$vz, '~', paste(edata$vg, collapse = '+'))
  if(length(edata$vx) > 0){
    eform <- paste(eform, '+', paste(edata$vx, collapse = '+'))
  }
  eform <- as.formula(eform)
  efit <- lm(eform, data = edata$data)
  efit

}
