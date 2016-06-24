
naive.method <- function(data.retro, data.expo){

  form <- create.formula(data.retro, data.expo)
  rform <- form$rform
  eform <- form$eform
  nform <- form$nform

  if(!is.null(nform)){
    afit <- lm(nform, data = data.expo$edat)
    data.expo$edat[, data.expo$expo.var] <- afit$residuals
    alp.a <- summary(afit)$coefficients[data.expo$add.covar.expo, 'Estimate']
    se.a <- summary(afit)$coefficients[data.expo$add.covar.expo, 'Std. Error']
    summary.a <- data.frame(Estimate = alp.a, SE = se.a, stringsAsFactors = FALSE)
    rownames(summary.a) <- paste0('alp.', data.expo$add.covar.expo)
    summary.a$z <- summary.a$Estimate / summary.a$SE
    summary.a$"Pr(>|z|)" <- pchisq(summary.a$z^2, df = 1, lower.tail = FALSE)
  }else{
    summary.a <- NULL
  }

  efit <- lm(eform, data = data.expo$edat)

  mu.y <- coef(efit)[data.expo$retro.var]

  alp.0 <- coef(efit)['(Intercept)']
  names(alp.0) <- 'alp.0'

  if(length(data.expo$overlap.covar) == 0){
    alp.o <- NA
    names(alp.o) <- 'alp.o'
  }else{
    alp.o <- coef(efit)[data.expo$overlap.covar]
    names(alp.o) <- paste0('alp.', data.expo$overlap.covar)
  }

  alp.g <- coef(efit)[data.expo$geno.var]
  names(alp.g) <- paste0('alp.', data.expo$geno.var)

  ne <- nrow(data.expo$edat)
  c <- log(mean(efit$residuals^2))
  names(c) <- 'c'

  data.retro$rdat$pred.expo <- as.matrix(data.retro$rdat[, data.retro$geno.var, drop = FALSE]) %*% alp.g

  rfit <- glm(rform, data = data.retro$rdat, family = 'binomial')

  bet.e <- coef(rfit)['pred.expo']
  names(bet.e) <- paste0('bet.', data.expo$expo.var)

  if(length(data.retro$add.covar.retro) > 0){
    bet.a <- coef(rfit)[data.retro$add.covar.retro]
    names(bet.a) <- paste0('bet.', data.retro$add.covar.retro)
  }else{
    bet.a <- NA
    names(bet.a) <- 'bet.a'
  }

  if(length(data.retro$overlap.covar) > 0){
    bet.o <- coef(rfit)[data.retro$overlap.covar] - bet.e * alp.o[paste0('alp.',data.retro$overlap.covar)]
    names(bet.o) <- paste0('bet.', data.retro$overlap.covar)
  }else{
    bet.o <- NA
    names(bet.o) <- 'bet.o'
  }

  b <- mu.y - exp(c) * bet.e
  names(b) <- 'b'

  if(!is.na(b)){
    a <- coef(rfit)['(Intercept)'] - (alp.0 + b) * bet.e - .5 * exp(c) * bet.e^2 # a = bet.0 + .5 * bet.u^2
  }else{
    a <- coef(rfit)['(Intercept)'] - alp.0 * bet.e - .5 * exp(c) * bet.e^2 # a = bet.0 + .5 * bet.u^2 + b * bet.z
  }
  names(a) <- 'a'

  list(naive.est = list(a = a, b = b, c = c, alp.0 = alp.0, alp.o = alp.o, alp.g = alp.g, bet.a = bet.a, bet.o = bet.o, bet.e = bet.e),
       summary.a = summary.a)

}

