

align.parameter <- function(fit, rm = FALSE){

  par.name <- NULL
  start <- NULL
  end <- NULL
  par <- NULL

  np <- length(fit)

  for(i in 1:np){
    if(rm && (names(fit)[i] == 'bet.z')){
      next
    }
    tmp <- fit[[names(fit)[i]]]
    if(all(!is.na(tmp))){
      par <- c(par, tmp)
      par.name <- c(par.name, names(fit)[i])
      rg <- range(which(names(par) %in% names(tmp)))
      start <- c(start, rg[1])
      end <- c(end, rg[2])
    }
  }

  par.pos <- data.frame(par.name, start, end, stringsAsFactors = FALSE)
  rownames(par.pos) <- par.name

  list(par = par, par.pos = par.pos)

}