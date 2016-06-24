

align.parameter <- function(naive.est){

  par.name <- NULL
  start <- NULL
  end <- NULL
  par <- NULL

  for(i in 1:9){
    tmp <- naive.est[[names(naive.est)[i]]]
    if(all(!is.na(tmp))){
      par <- c(par, tmp)
      par.name <- c(par.name, names(naive.est)[i])
      rg <- range(which(names(par) %in% names(tmp)))
      start <- c(start, rg[1])
      end <- c(end, rg[2])
    }
  }

  par.pos <- data.frame(par.name, start, end, stringsAsFactors = FALSE)
  rownames(par.pos) <- par.name

  list(par = par, par.pos = par.pos)

}