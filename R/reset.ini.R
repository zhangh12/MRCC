
reset.ini <- function(par, nseed){

  message('reset initial guess')
  if(nseed == 0){
    par[] <- 0
  }else{
    par[] <- par + runif(length(par), -.2, .2)
  }

  par

}