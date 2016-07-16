
naive.method <- function(rdata, edata){

  form <- create.formula(rdata, edata)
  rform <- form$rform
  eform <- form$eform
  nform <- form$nform

  if(!is.null(nform)){
    hfit <- lm(nform, data = edata$data)
    edata$data[, edata$vz] <- hfit$residuals
    alp.h <- summary(hfit)$coefficients[edata$vh, 'Estimate']
    se.h <- summary(hfit)$coefficients[edata$vh, 'Std. Error']
    summary.h <- data.frame(Estimate = alp.h, SE = se.h, SE0 = se.h, stringsAsFactors = FALSE)
    rownames(summary.h) <- paste0('alp.', edata$vh)
    summary.h$z <- summary.h$Estimate / summary.h$SE
    summary.h$"Pr(>|z|)" <- pchisq(summary.h$z^2, df = 1, lower.tail = FALSE)
  }else{
    summary.h <- NULL
  }

  efit <- lm(eform, data = edata$data)

  mu.d <- coef(efit)[edata$vd]

  alp.0 <- coef(efit)['(Intercept)']
  names(alp.0) <- 'alp.0'

  if(length(edata$vx) == 0){
    alp.x <- NA
    names(alp.x) <- 'alp.x'
  }else{
    alp.x <- coef(efit)[edata$vx]
    names(alp.x) <- paste0('alp.', edata$vx)
  }

  alp.g <- coef(efit)[edata$vg]
  names(alp.g) <- paste0('alp.', edata$vg)

  ne <- nrow(edata$data)
  c <- log(mean(efit$residuals^2))
  names(c) <- 'c'

  rdata$data$pred.expo <- as.matrix(rdata$data[, rdata$vg, drop = FALSE]) %*% alp.g

  rfit <- glm(rform, data = rdata$data, family = 'binomial')

  bet.z <- coef(rfit)['pred.expo']
  names(bet.z) <- paste0('bet.', edata$vz)

  if(length(rdata$vy) > 0){
    bet.y <- coef(rfit)[rdata$vy]
    names(bet.y) <- paste0('bet.', rdata$vy)
  }else{
    bet.y <- NA
    names(bet.y) <- 'bet.y'
  }

  if(length(rdata$vx) > 0){
    bet.x <- coef(rfit)[rdata$vx] - bet.z * alp.x[paste0('alp.',rdata$vx)]
    names(bet.x) <- paste0('bet.', rdata$vx)
  }else{
    bet.x <- NA
    names(bet.x) <- 'bet.x'
  }

  b <- mu.d - exp(c) * bet.z
  names(b) <- 'b'

  if(!is.na(b)){
    a <- coef(rfit)['(Intercept)'] - (alp.0 + b) * bet.z - .5 * exp(c) * bet.z^2 # a = bet.0 + .5 * bet.u^2
  }else{
    a <- coef(rfit)['(Intercept)'] - alp.0 * bet.z - .5 * exp(c) * bet.z^2 # a = bet.0 + .5 * bet.u^2 + b * bet.z
  }
  names(a) <- 'a'

  list(naive.est = list(a = a, b = b, c = c, alp.0 = alp.0, alp.x = alp.x, alp.g = alp.g, bet.y = bet.y, bet.x = bet.x, bet.z = bet.z),
       summary.h = summary.h, edata = edata)

}

