
#' Generate trancated random normal variables where values outside range are assigned to range limit
#'
#' Note that this function is different from the one in the package truncnorm 
#'
#' Updated 2018-08-14

#' @param n number of variables to be generated
#' @param a lower limit of range
#' @param b upper limit of range
#' @param meana mean of underlying normal distribution prior to truncation
#' @param sda standard deviation of underlying normal distribution prior to truncation
#' @keywords misc
#' @export 
#' @examples
#' truncnorm()

# to do - GENERAL TESTING 


# for generating truncated normal random variables
altrtruncnorm <- function(n,a=0,b=1,meana=0,sda=1){
  j <- rnorm(n,mean=meana,sd=sda)
j[j < a] <- a
j[j > b] <- b
j
}
