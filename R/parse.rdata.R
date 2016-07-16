
parse.rdata <- function(rformula, rdata){

  if(!("Formula" %in% class(rformula))){
    if("formula" %in% class(rformula)){
      rformula<-Formula(rformula)
    }else{
      stop("rformula should be of class \"formula\"")
    }
  }

  if(!('id' %in% colnames(rdata))){
    stop('Cannot find a column \'id\' in rdata')
  }

  if(any(duplicated(rdata$id))){
    stop('ID in rdata is not unique')
  }

  rownames(rdata) <- rdata$id

  mf <- model.frame(rformula, na.action = na.pass, data = rdata, rhs=1:2, lhs=1, drop=FALSE)
  covar <- model.matrix(rformula, mf, rhs=1, drop=F)[, -1, drop = FALSE]
  geno <- model.matrix(rformula, mf, rhs=2, drop=F)[, -1, drop = FALSE]

  retro <- model.part(rformula, mf, lhs=1, drop=FALSE)

  retro.var <- colnames(retro)
  covar.var <- colnames(covar)
  geno.var <- colnames(geno)
  data <- data.frame(retro, covar, geno, stringsAsFactors = FALSE)

  list(data = data, retro.var = retro.var, covar.var = covar.var, geno.var = geno.var)

}
