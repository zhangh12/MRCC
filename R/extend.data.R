
extend.data <- function(data.retro, data.expo){

  overlap.covar <- intersect(data.expo$covar.var, data.retro$covar.var)
  if(is.null(overlap.covar)){
    overlap.covar <- character()
  }
  data.expo$overlap.covar <- overlap.covar
  data.retro$overlap.covar <- overlap.covar

  add.covar.expo <- setdiff(data.expo$covar.var, data.retro$covar.var)
  if(is.null(add.covar.expo)){
    add.covar.expo <- character()
  }
  data.expo$add.covar.expo <- add.covar.expo

  for(v in add.covar.expo){
    data.retro$rdat[, v] <- 0
  }

  add.covar.retro <- setdiff(data.retro$covar.var, data.expo$covar.var)
  if(is.null(add.covar.retro)){
    add.covar.retro <- character()
  }
  data.retro$add.covar.retro <- add.covar.retro

  for(v in add.covar.retro){
    data.expo$edat[, v] <- 0
  }

  if(setequal(data.retro$geno.var, data.expo$geno.var)){
    data.retro$geno.var <- data.expo$geno.var
  }else{
    stop('The same set of SNPs should be measured in case-control study and exposure data. ')
  }

  list(data.retro = data.retro, data.expo = data.expo)


}
