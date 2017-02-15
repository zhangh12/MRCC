
LRT.conf.int <- function(rdata, edata, logL.obs.alt, par, se, c.lrt, level, fig){

  crit <- c.lrt * qchisq(level, df = 1)
  par.pos <- align.parameter(find.tsls(rdata, edata))$par.pos

  name.bet.z <- paste0('bet.', edata$vz)
  se <- se[name.bet.z]
  b0 <- par[name.bet.z]
  names(b0) <- NULL
  names(se) <- NULL
  rm(par)
  bet <- seq(b0 - 10 * qnorm((1+level)/2) * se, b0 + 10 * qnorm((1+level)/2) * se, by = .05)
  bet <- c(bet, 0)
  bet[abs(bet) < 1e-3] <- 0
  bet <- sort(unique(bet))

  lr.stat <- NULL

  for(b in bet){
    if(b == 0){
      par.mle <- align.parameter(find.null.tsls(rdata, edata))$par
    }else{
      par.mle <-  try(find.LRT.mle(rdata, edata, b), silent = TRUE)
      if('try-error' %in% class(par.mle)){
        lr.stat <- c(lr.stat, NA)
        next
      }
      par.mle <- c(par.mle, b)
      names(par.mle)[par.pos[name.bet.z, 'start']] <- name.bet.z
    }

    lr.stat <- c(lr.stat, 2 * (logL.obs.alt - logL(par.mle, rdata, edata, par.pos)))
    if(fig){
      plot(bet[1:length(lr.stat)],lr.stat,pch=20,cex=.2,type='b',ylim=range(c(crit,lr.stat), na.rm = TRUE))
      abline(h=crit)
      abline(v=b0)
    }
  }


  bet <- bet[!is.na(lr.stat)]
  lr.stat <- na.omit(lr.stat)

  delta0 <- lr.stat - crit
  delta1 <- head(delta0, -1)
  delta2 <- delta0[-1]

  # approx points for left CI
  left1 <- NULL
  left2 <- NULL

  if(delta1[1] < 0){
    left1 <- -Inf
    left2 <- -Inf
  }

  id <- which(delta1 >= 0 & delta2 <= 0)
  if(length(id) > 0){
    left1 <- c(left1, bet[id])
    left2 <- c(left2, bet[id+1])
  }

  # approx points for right CI
  right1 <- NULL
  right2 <- NULL

  id <- which(delta1 <= 0 & delta2 >= 0)
  if(length(id) > 0){
    right1 <- c(right1, bet[id])
    right2 <- c(right2, bet[id+1])
  }

  if(tail(delta0, 1) < 0){
    right1 <- c(right1, Inf)
    right2 <- c(right2, Inf)
  }

  if(length(left1) != length(right1)){
    stop('debug')
  }

  nci <- length(left1)

  # search for left end point of CI
  lci <- NULL
  for(i in 1:nci){
    b1 <- left1[i]
    b2 <- left2[i]

    if(is.infinite(b1)){
      lci <- c(lci, b1)
      next
    }

    bs <- NULL
    while(TRUE){
      b <- (b1 + b2)/2
      bs <- c(bs, b)
      cond <- ((length(bs) >2 && sd(tail(bs, 2)) == 0) || (length(bs) > 20 && sd(tail(bs, 10))  < 1e-3))
      if(cond){
        lci <- c(lci, b)
        break
      }

      if(abs(b) < 1e-3){
        par.mle <- align.parameter(find.null.tsls(rdata, edata))$par
      }else{
        par.mle <-  try(find.LRT.mle(rdata, edata, b), silent = TRUE)
        if('try-error' %in% class(par.mle)){
          stop('error in LRT.conf.int')
        }
        par.mle <- c(par.mle, b)
        names(par.mle)[par.pos[name.bet.z, 'start']] <- name.bet.z
      }

      stat <- 2 * (logL.obs.alt - logL(par.mle, rdata, edata, par.pos))
      if(stat > crit + 1e-3){
        b1 <- b
        if(fig){
          plot(bet,lr.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      if(stat < crit - 1e-3){
        b0 <- b
        if(fig){
          plot(bet,lr.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      lci <- c(lci, b)
      if(fig){
        plot(bet,lr.stat,pch=20,cex=.2,type='b')
        abline(v = b,col='red')
        abline(h=crit,col='blue')
      }
      break
    }
  }

  # search for right end point of CI
  rci <- NULL
  for(i in 1:nci){
    b1 <- right1[i]
    b2 <- right2[i]

    if(is.infinite(b1)){
      rci <- c(rci, b1)
      next
    }

    bs <- NULL
    while(TRUE){
      b <- (b1 + b2)/2
      bs <- c(bs, b)
      cond <- ((length(bs) >2 && sd(tail(bs, 2)) == 0) || (length(bs) > 20 && sd(tail(bs, 10))  < 1e-3))
      if(cond){
        rci <- c(rci, b)
        break
      }

      if(abs(b) < 1e-3){
        par.mle <- align.parameter(find.null.tsls(rdata, edata))$par
      }else{
        par.mle <-  try(find.LRT.mle(rdata, edata, b), silent = TRUE)
        if('try-error' %in% class(par.mle)){
          stop('error in LRT.conf.int')
        }
        par.mle <- c(par.mle, b)
        names(par.mle)[par.pos[name.bet.z, 'start']] <- name.bet.z
      }

      stat <- 2 * (logL.obs.alt - logL(par.mle, rdata, edata, par.pos))
      if(stat > crit + 1e-3){
        b2 <- b
        if(fig){
          plot(bet,lr.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      if(stat < crit - 1e-3){
        b1 <- b
        if(fig){
          plot(bet,lr.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      rci <- c(rci, b)
      if(fig){
        plot(bet,lr.stat,pch=20,cex=.2,type='b')
        abline(v = b,col='red')
        abline(h=crit,col='blue')
      }
      break
    }
  }

  data.frame(LCL = lci, RCL = rci)

}