

generate.null.data <- function(rdata, edata, mdata, efit){

  er0 <- efit$residuals

  id3 <- intersect(rownames(rdata$data), rownames(edata$data))
  id1 <- setdiff(rownames(rdata$data), id3)
  id2 <- setdiff(rownames(edata$data), id3)

  id13 <- c(id1, id3)
  id23 <- c(id2, id3)

  nid13 <- sample(id13)

  vx <- rdata$vx
  vy <- rdata$vy
  vd <- rdata$vd
  vg <- rdata$vg
  vz <- edata$vz

  mdata[id13, vg] <- mdata[nid13, vg, drop = FALSE]

  n23 <- length(id23)
  pred <- predict.lm(efit, mdata[id23, c(vx, vg)])
  err <- sample(er0, n23, TRUE)
  mdata[id23, vz] <- pred + err

  if(!is.null(vx)){
    rx <- mdata[id13, vx, drop = FALSE]
    ex <- mdata[id23, vx, drop = FALSE]
  }else{
    rx <- NULL
    ex <- NULL
  }

  if(!is.null(vy)){
    ry <- mdata[id13, vy, drop = FALSE]
  }else{
    ry <- NULL
  }

  rg <- mdata[id13, vg, drop = FALSE]
  eg <- mdata[id23, vg, drop = FALSE]

  rd <- mdata[id13, vd, drop = FALSE]
  ez <- mdata[id23, vz, drop = FALSE]

  rdat <- mdata[id13, c(vd, vg, vx, vy), drop = FALSE]
  edat <- mdata[id23, c(vz, vg, vx), drop = FALSE]

  rdata <- list(data = rdat,
                vx = vx, vy = vy, vd = vd, vg = vg,
                rx = rx, ry = ry, rd = rd, rg = rg)
  edata <- list(data = edat,
                vx = vx, vg = vg, vz = vz,
                ex = ex, eg = eg, ez = ez)

  list(rdata = rdata, edata = edata)

}
