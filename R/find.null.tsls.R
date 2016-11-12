
find.null.tsls <- function(rdata, edata){

  form <- create.null.formula(rdata, edata)
  rform <- form$rform
  eform <- form$eform

  efit <- lm(eform, data = edata$data)

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
  c <- mean(efit$residuals^2)
  names(c) <- 'c'

  rfit <- glm(rform, data = rdata$data, family = 'binomial')

  if(length(rdata$vy) > 0){
    bet.y <- coef(rfit)[rdata$vy]
    names(bet.y) <- paste0('bet.', rdata$vy)
  }else{
    bet.y <- NA
    names(bet.y) <- 'bet.y'
  }

  if(length(rdata$vx) > 0){
    bet.x <- coef(rfit)[rdata$vx]
    names(bet.x) <- paste0('bet.', rdata$vx)
  }else{
    bet.x <- NA
    names(bet.x) <- 'bet.x'
  }

  n1 <- sum(rdata$data[, rdata$vd])
  n0 <- nrow(rdata$data) - n1
  a <- coef(rfit)['(Intercept)'] - log(n1/n0)
  names(a) <- 'a'

  bet.z <- 0
  names(bet.z) <- paste0('bet.', edata$vz)

  list(c = c, alp.0 = alp.0, alp.x = alp.x, alp.g = alp.g, a = a, bet.x = bet.x, bet.y = bet.y, bet.z = bet.z)

}
