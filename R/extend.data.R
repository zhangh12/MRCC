
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
  ed <- as.matrix(edata[, vd, drop = FALSE])

  rdat <- rdata[, c(vd, vg, vx, vy)]
  edat <- edata[, c(vd, vg, vx, vh, vz)]

  rdata <- list(data = rdat,
                vx = vx, vy = vy, vd = vd, vg = vg,
                rx = rx, ry = ry, rd = rd, rg = rg)
  edata <- list(data = edat,
                vx = vx, vg = vg, vz = vz, vd = vd, vh = vh,
                ex = ex, eg = eg, ez = ez, ed = ed, eh = eh)

  id1 <- setdiff(rownames(rdat), rownames(edat))
  id2 <- setdiff(rownames(edat), rownames(rdat))
  id3 <- intersect(rownames(rdat), rownames(edat))

  id10 <- id1[which(rdat[id1, vd] == 0)]
  id11 <- id1[which(rdat[id1, vd] == 1)]

  id20 <- id2[which(edat[id2, vd] == 0)]
  id21 <- id2[which(edat[id2, vd] == 1)]

  id30 <- id3[which(rdat[id3, vd] == 0)]
  id31 <- id3[which(rdat[id3, vd] == 1)]

  omega <- list(id1 = id1, id10 = id10, id11 = id11,
                id2 = id2, id20 = id20, id21 = id21,
                id3 = id3, id30 = id30, id31 = id31)


  list(rdata = rdata, edata = edata, omega = omega)

}
