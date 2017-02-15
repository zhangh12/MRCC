
LMT.conf.int <- function(rdata, edata, par, se, c.lmt, level, fig){

  crit <- qchisq(level, df = 1)
  par.pos <- align.parameter(find.null.tsls(rdata, edata))$par.pos

  name.bet.z <- paste0('bet.', edata$vz)
  se <- se[name.bet.z]
  b0 <- par[name.bet.z]
  names(b0) <- NULL
  names(se) <- NULL
  rm(par)
  bet <- seq(b0 - 20 * qnorm((1+level)/2) * se, b0 + 20 * qnorm((1+level)/2) * se, by = .1)
  bet <- c(bet, 0)
  bet[abs(bet) < 1e-3] <- 0
  bet <- sort(unique(bet))

  lm.stat <- NULL
  left1 <- -Inf
  left2 <- Inf
  right1 <- -Inf
  right2 <- Inf

  for(b in bet){
    if(abs(b) < 1e-3){
      par.mle <- align.parameter(find.null.tsls(rdata, edata))$par
    }else{
      par.mle <-  try(find.LRT.mle(rdata, edata, b), silent = TRUE)
      if('try-error' %in% class(par.mle)){
        lm.stat <- c(lm.stat, NA)
        next
      }
      par.mle <- c(par.mle, b)
      names(par.mle)[par.pos[name.bet.z, 'start']] <- name.bet.z
    }

    J0 <- -hessian0(par.mle, rdata, edata, par.pos)
    I0 <- fisher.info0(par.mle, rdata, edata, par.pos)
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

    tmp.b <- bet[!is.na(lm.stat)]
    tmp <- lm.stat[!is.na(lm.stat)]
    len <- length(tmp)
    if(len > 2){
      if(tmp[len] <= crit && tmp[len-1] >= crit){
        left1 <- tmp.b[len-1]
        right1 <- tmp.b[len]
        next
      }

      if(tmp[len] >= crit && tmp[len-1] <= crit){
        left2 <- tmp.b[len-1]
        right2 <- tmp.b[len]
        #break
      }
    }
  }

  bet <- bet[1:length(lm.stat)]

  if(is.infinite(left1)){
    lci <- -Inf
  }else{
    lci <- NULL
  }

  if(is.infinite(left2)){
    rci <- Inf
  }else{
    rci <- NULL
  }

  if(any(is.infinite(c(lci, rci)))){
    return(c(LCL = NA, RCL = NA))
  }

  # search for left end point of CI
  b00 <- NULL
  b0 <- right1
  b1 <- left1
  if(is.null(lci)){
    while(TRUE){
      b <- (b0 + b1)/2
      b00 <- c(b00, b)
      if(length(b00) > 20 && sd(tail(b00, 10)) < 1e-3){
        lci <- b
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
        names(par.mle)[par.pos[name.bet.z, 'start']] <- name.bet.z
      }

      J0 <- -hessian0(par.mle, rdata, edata, par.pos)
      I0 <- fisher.info0(par.mle, rdata, edata, par.pos)
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
        b0 <- b
        if(fig){
          plot(bet,lm.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      lci <- b
      break
    }
  }

  # search for right end point of CI
  b00 <- NULL
  b0 <- left2
  b2 <- right2
  if(is.null(rci)){
    while(TRUE){
      b <- (b0 + b2)/2
      b00 <- c(b00, b)
      if(length(b00) > 20 && sd(tail(b00, 10)) < 1e-3){
        rci <- b
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
        names(par.mle)[par.pos[name.bet.z, 'start']] <- name.bet.z
      }

      J0 <- -hessian0(par.mle, rdata, edata, par.pos)
      I0 <- fisher.info0(par.mle, rdata, edata, par.pos)
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
        b0 <- b
        if(fig){
          plot(bet,lm.stat,pch=20,cex=.2,type='b')
          abline(v = b,col='red')
          abline(h=crit,col='blue')
        }
        next
      }

      rci <- b
      break
    }
  }

  c(LCL = lci, RCL = rci)

}
