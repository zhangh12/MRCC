
extend.data <- function(rdata, edata){

  vx <- intersect(edata$covar.var, rdata$covar.var)
  if(is.null(vx)){
    vx <- character()
  }
  edata$vx <- vx
  rdata$vx <- vx

  vh <- setdiff(edata$covar.var, rdata$covar.var)
  if(is.null(vh)){
    vh <- character()
  }
  edata$vh <- vh

  for(v in vh){
    rdata$data[, v] <- 0
  }

  vy <- setdiff(rdata$covar.var, edata$covar.var)
  if(is.null(vy)){
    vy <- character()
  }
  rdata$vy <- vy

  for(v in vy){
    edata$data[, v] <- 0
  }

  if(setequal(rdata$geno.var, edata$geno.var)){
    rdata$vg <- rdata$geno.var
    edata$vg <- rdata$geno.var
  }else{
    stop('The same set of SNPs should be measured in case-control study and exposure data. ')
  }

  rdata$vd <- rdata$retro.var
  rdata$vg <- rdata$geno.var
  edata$vd <- edata$retro.var
  edata$vz <- edata$expo.var
  edata$vg <- edata$geno.var

  rdata <- list(data = rdata$data, vx = rdata$vx, vy = rdata$vy, vd = rdata$vd, vg = rdata$vg)
  edata <- list(data = edata$data, vx = edata$vx, vh = edata$vh, vd = edata$vd, vg = edata$vg, vz = edata$vz)

  list(rdata = rdata, edata = edata)


}
