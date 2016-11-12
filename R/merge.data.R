
merge.data <- function(rdata, edata){

  vx <- rdata$vx
  vy <- rdata$vy
  vd <- rdata$vd
  vg <- rdata$vg
  vz <- edata$vz

  rdat <- rdata$data
  edat <- edata$data
  rdat$id <- rownames(rdat)
  edat$id <- rownames(edat)

  mdata <- merge(rdat[, c(vd, vx, vy, vg, 'id')], edat[, c(vx, vg, vz, 'id')], by = c('id', vx, vg), all = TRUE)
  rownames(mdata) <- mdata$id
  mdata$id <- NULL

  mdata

}