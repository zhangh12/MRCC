
extend.data <- function(rdata, edata){

  vx <- intersect(edata$covar.var, rdata$covar.var)
  if(length(vx) == 0){
    vx <- NULL
  }

  vh <- setdiff(edata$covar.var, rdata$covar.var)
  if(length(vh) == 0){
    vh <- NULL
  }

  vy <- setdiff(rdata$covar.var, edata$covar.var)
  if(length(vy) == 0){
    vy <- NULL
  }

  if(setequal(rdata$geno.var, edata$geno.var)){
    vg <- rdata$geno.var
  }else{
    stop('The same set of SNPs should be measured in case-control study and exposure data. ')
  }

  vd <- rdata$retro.var
  vz <- edata$expo.var

  rdata <- rdata$data
  edata <- edata$data

  if(!is.null(vx)){
    rx <- as.matrix(rdata[, vx, drop = FALSE])
    ex <- as.matrix(edata[, vx, drop = FALSE])
  }else{
    rx <- NULL
    ex <- NULL
  }

  if(!is.null(vy)){
    ry <- as.matrix(rdata[, vy, drop = FALSE])
  }else{
    ry <- NULL
  }

  if(!is.null(vh)){
    eh <- as.matrix(edata[, vh, drop = FALSE])
  }else{
    eh <- NULL
  }

  rd <- as.matrix(rdata[, vd, drop = FALSE])

  rg <- as.matrix(rdata[, vg, drop = FALSE])
  eg <- as.matrix(edata[, vg, drop = FALSE])
  ez <- as.matrix(edata[, vz, drop = FALSE])

  rdat <- rdata[, c(vd, vg, vx, vy)]
  edat <- edata[, c(vz, vg, vx, vh)]

  if(!is.null(vh)){
    hform <- paste(vz, '~', paste(vh, collapse = '+'))
    hform <- as.formula(hform)
    hfit <- lm(hform, data = edat)
    ez[] <- hfit$residuals
    edat[, vz] <- ez
    alp.h <- summary(hfit)$coefficients[vh, 'Estimate']
    se.h <- summary(hfit)$coefficients[vh, 'Std. Error']
    summary.h <- data.frame(Estimate = alp.h, SE = se.h, SE0 = se.h, stringsAsFactors = FALSE)
    rownames(summary.h) <- paste0('alp.', vh)
    summary.h$z <- summary.h$Estimate / summary.h$SE
    summary.h$"Pr(>|z|)" <- pchisq(summary.h$z^2, df = 1, lower.tail = FALSE)
    edat[, vh] <- NULL
    vh <- NULL
  }else{
    summary.h <- NULL
  }

  rdata <- list(data = rdat,
                vx = vx, vy = vy, vd = vd, vg = vg,
                rx = rx, ry = ry, rd = rd, rg = rg)
  edata <- list(data = edat,
                vx = vx, vg = vg, vz = vz,
                ex = ex, eg = eg, ez = ez)

  list(rdata = rdata, edata = edata, summary.h = summary.h)

}
