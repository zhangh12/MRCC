
null.adjustment <- function(rdata, edata, nb = 500){

  efit <- expo.model(edata)
  mdata <- merge.data(rdata, edata)

  lrt.b <- NULL
  wald.b <- NULL
  for(k in 1:max(1000, 2*nb)){
    set.seed(k)
    null.data <- generate.null.data(rdata, edata, mdata, efit)
    t1 <- try(wald <- wald.test(null.data$rdata, null.data$edata, 1.0, 0.95, TRUE), silent = TRUE)
    if('try-error' %in% class(t1)){
      next
    }

    logL.b.alt <- wald$max.logL
    wald.b <- c(wald.b, wald$wald.stat)
    logL.b.null <- logL.null(null.data$rdata, null.data$edata)
    lrt.b <- c(lrt.b, -2 * (logL.b.null - logL.b.alt))

    if(FALSE){
      rdat <- null.data$rdata$data
      edat <- null.data$edata$data
      rd <- rdat[, null.data$rdata$vd]
      rg <- rdat[, null.data$rdata$vg]
      ez <- edat[, null.data$edata$vz]
      eg <- edat[, null.data$edata$vg]
      save(list=c('rd', 'rg', 'ez', 'eg'), file = '/Users/zhangh12/test.rda')
    }

    #print(k)
    if(length(lrt.b) == nb){
      break
    }
  }

  c.lrt <- mean(lrt.b)
  c.wald <- mean(wald.b)

  list(c.lrt = c.lrt, c.wald = c.wald)

}

