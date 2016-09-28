
LRT <- function(rdata, edata, fn){

  id3 <- intersect(rownames(rdata$data), rownames(edata$data))

  id10 <- setdiff(rownames(rdata$data[rdata$data[, rdata$vd] == 0, , drop = FALSE]), id3)
  id11 <- setdiff(rownames(rdata$data[rdata$data[, rdata$vd] == 1, , drop = FALSE]), id3)

  n0 <- length(c(id10, id3))
  n1 <- length(id11)
  ne <- nrow(edata$data)

  rform <- paste0(rdata$vd, ' ~ ')
  vxy <- c(rdata$vx, rdata$vy)
  if(length(vxy) > 0){
    rform <- paste0(rform, paste0(vxy, collapse = ' + '))
  }else{
    rform <- paste0(rform, '1')
  }

  eform <- paste0(edata$vz, ' ~ ')
  vx <- edata$vx
  if(length(vx) > 0){
    eform <- paste0(eform, paste0(vx, collapse = ' + '))
  }else{
    eform <- paste0(eform, '1')
  }

  rfit <- glm(rform, data = rdata$data, family = 'binomial')
  rhat <- rfit$fitted.values
  rlin <- log(rhat/(1-rhat))

  rd <- rdata$data[, rdata$vd, drop = TRUE]

  rfn <- -sum(rd * rlin) + sum(log(1 + exp(rlin))) + n1*log(n1)+n0*log(n0)

  efit <- lm(eform, data = edata$data)
  eres <- efit$residuals
  c <- log(mean(eres^2))
  efn <- ne/2 * log(2 * pi) + ne/2 * c + ne/2

  fn0 <- rfn + efn

  lrt <- 2 * (fn0 - fn)

  pv <- pchisq(lrt, df = 1, lower.tail = FALSE)

  pv

}
