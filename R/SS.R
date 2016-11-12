

SS <- function(rd,rg,ez,eg){

  b.x <- NULL
  se.x <- NULL
  b.y <- NULL
  se.y <- NULL
  np <- ncol(rg)
  for(i in 1:np){
    efit <- lm(paste0('z~g',i), data.frame(z=ez,eg))
    rfit <- glm(paste0('y~g',i),data.frame(y=rd,rg),family='binomial')
    b.x <- c(b.x, coef(efit)[2])
    b.y <- c(b.y, coef(rfit)[2])
    se.x <- c(se.x, summary(efit)$coefficients[2, 'Std. Error'])
    se.y <- c(se.y, summary(rfit)$coefficients[2, 'Std. Error'])
  }

  bet.b <- sum(b.x*b.y/se.y^2)/sum(b.x^2/se.y^2) # eq (4) Burgess et al. (2015) Stat Med
  se.b <- sqrt(1/sum(b.x^2/se.y^2)) # eq (5) Burgess et al. (2015) Stat Med


  # Mokry et al. (2015) PLoS Med, exactly the same as the inverse-variance weighted in Burgess et al. (2015) Stat Med
  # w <- 1/(se.y/b.x)^2
  # bet.meta <- sum(w*b.y/b.x)/sum(w) # eq (1) Mokry et al. (2015) PLoS Med
  # se.meta <- sqrt(1/sum(w)) # eq (2)  Mokry et al. (2015) PLoS Med

  #causal effect, yang's method
  se <- sqrt(b.y^2/b.x^2* (se.x^2/b.x^2 + se.y^2/b.y^2))
  bet.d <- sum(b.y/b.x/se^2)/sum(1/se^2)
  se.d <- sqrt(1/sum(1/se^2))

  c(bet.b =bet.b,se.b=se.b,bet.d=bet.d,se.d=se.d)

}
