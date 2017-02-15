

align.parameter <- function(fit, null = FALSE){

  par.name <- NULL
  start <- NULL
  end <- NULL
  par <- NULL

  np <- ifelse(null, 7, 8)

  for(i in 1:np){
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