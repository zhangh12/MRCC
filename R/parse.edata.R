
parse.edata <- function(eformula, edata){

  if(!("Formula" %in% class(eformula))){
    if("formula" %in% class(eformula)){
      eformula<-Formula(eformula)
    }else{
      stop("rformula should be of class \"formula\"")
    }
  }

  if(!('id' %in% colnames(edata))){
    stop('Cannot find a column \'id\' in edata')
  }

  if(any(duplicated(edata$id))){
    stop('ID in edata is not unique')
  }

  rownames(edata) <- edata$id

  mf <- model.frame(eformula, na.action = na.pass, data = edata, rhs=1:2, lhs=1, drop=FALSE)
  covar <- model.matrix(eformula, mf, rhs=1, drop=F)[, -1, drop = FALSE]
  geno <- model.matrix(eformula, mf, rhs=2, drop=F)[, -1, drop = FALSE]

  expo <- model.part(eformula, mf, lhs=1, drop=FALSE)

  expo.var <- colnames(expo)
  covar.var <- colnames(covar)
  geno.var <- colnames(geno)

  data <- data.frame(expo, covar, geno, stringsAsFactors = FALSE)

  list(data = data, expo.var = expo.var, covar.var = covar.var, geno.var = geno.var)

}
