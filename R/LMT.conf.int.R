
LMT.conf.int <- function(rdata, edata, omega, par, se, level, fig){

  crit <- qchisq(level, df = 1)
  par.pos <- align.parameter(find.null.tsls(rdata, edata))$par.pos

  name.bet.z <- paste0('bet.', edata$vz)
  se <- se[name.bet.z]
  b0 <- par[name.bet.z]
  names(b0) <- NULL
  names(se) <- NULL
  rm(par)
  bet <- seq(b0 - 5 * qnorm((1+level)/2) * se, b0 + 5 * qnorm((1+level)/2) * se, by = .05)
  bet <- c(bet, 0)
  bet[abs(bet) < 1e-3] <- 0
  bet <- sort(unique(bet))

  lm.stat <- NULL

  for(b in bet){
    if(b == 0){
      par.mle <- align.parameter(find.null.tsls(rdata, edata))$par
    }else{
      par.mle <-  try(find.LRT.mle(rdata, edata, b), silent = TRUE)
      if('try-error' %in% class(par.mle)){
        lm.stat <- c(lm.stat, NA)
        next
      }
      par.mle <- c(par.mle, b)
      names(par.mle)[par.pos['bet.z', 'start']] <- name.bet.z
    }

    J0 <- -hessian0(par.mle, rdata, edata, omega, par.pos)
    I0 <- fisher.info0(par.mle, rdata, edata, omega, par.pos)
    sc <- score(par.mle, rdata, edata, par.pos)[name.bet.z]

    id.bet <- which(colnames(J0) == name.bet.z)
    v <- cbind(-J0[id.bet,-id.bet] %*% solve(J0[-id.bet, -id.bet]), 1)
    colnames(v)[ncol(v)] <- name.bet.z
    v <- v[, colnames(I0), drop = FALSE]

    lm.stat <- c(lm.stat, sc^2/(v %*% I0 %*% t(v)))
    if(fig){
      plot(bet[1:length(lm.stat)],lm.stat,pch=20,cex=.2,type='b',ylim=range(c(crit,lm.stat), na.rm = TRUE))
      abline(h=crit)
      abline(v=b0)
    }

  }

  bet <- bet[!is.na(lm.stat)]
  lm.stat <- na.omit(lm.stat)

  delta0 <- lm.stat - crit
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
          stop('error in LMT.conf.int')
        }
        par.mle <- c(par.mle, b)
        names(par.mle)[par.pos['bet.z', 'start']] <- name.bet.z
      }

      J0 <- -hessian0(par.mle, rdata, edata, omega, par.pos)
      I0 <- fisher.info0(par.mle, rdata, edata, omega, par.pos)
      sc <- score(par.mle, rdata, edata, par.pos)[name.bet.z]
      id.bet <- which(colnames(J0) == name.bet.z)
      v <- cbind(-J0[id.bet,-id.bet] %*% solve(J0[-id.bet, -id.bet]), 1)
      colnames(v)[ncol(v)] <- name.bet.z
      v <- v[, colnames(I0), drop = FALSE]

      stat <- sc^2/(v %*% I0 %*% t(v))
      if(stat > crit + 1e-3){
        b1 <- b
        if(fig){
          plot(bet,lm.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      if(stat < crit - 1e-3){
        b2 <- b
        if(fig){
          plot(bet,lm.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      lci <- c(lci, b)
      if(fig){
        plot(bet,lm.stat,pch=20,cex=.2,type='b')
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
          stop('error in LMT.conf.int')
        }
        par.mle <- c(par.mle, b)
        names(par.mle)[par.pos['bet.z', 'start']] <- name.bet.z
      }

      J0 <- -hessian0(par.mle, rdata, edata, omega, par.pos)
      I0 <- fisher.info0(par.mle, rdata, edata, omega, par.pos)
      sc <- score(par.mle, rdata, edata, par.pos)[name.bet.z]
      id.bet <- which(colnames(J0) == name.bet.z)
      v <- cbind(-J0[id.bet,-id.bet] %*% solve(J0[-id.bet, -id.bet]), 1)
      colnames(v)[ncol(v)] <- name.bet.z
      v <- v[, colnames(I0), drop = FALSE]

      stat <- sc^2/(v %*% I0 %*% t(v))
      if(stat > crit + 1e-3){
        b2 <- b
        if(fig){
          plot(bet,lm.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      if(stat < crit - 1e-3){
        b1 <- b
        if(fig){
          plot(bet,lm.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      rci <- c(rci, b)
      if(fig){
        plot(bet,lm.stat,pch=20,cex=.2,type='b')
        abline(v = b,col='red')
        abline(h=crit,col='blue')
      }
      break
    }
  }

  #id <- which(b0 >= lci & b0 <= rci)
  data.frame(LCL = lci, RCL = rci)

}
