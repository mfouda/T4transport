## Auxiliary Computations
## (1) create_common_probability



# (1) create_common_probability -------------------------------------------
#     Given a length-K list of weights, create a common matrix where
#     each column is an atomic weights.
#' @keywords internal
create_common_probability <- function(weights){
  K = length(weights)
  mysize = unlist(lapply(weights, length))
  allsize = sum(mysize)
  
  output = array(0,c(allsize, K))
  for (k in 1:K){
    if (k==1){
      output[,k] = c(weights[[k]], rep(0,sum(mysize[(k+1):K])))
    } else if (k==K){
      output[,k] = c(rep(0,sum(mysize[1:(K-1)])), weights[[k]])
    } else {
      output[,k] = c(rep(0,sum(mysize[1:(k-1)])), weights[[k]], rep(0,sum(mysize[(k+1):K])))
    }
  }
  return(output)
}
