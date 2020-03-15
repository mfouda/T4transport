#' Barycenter of Empirical Measures with Fixed Support via Proximal Algorithm
#' 
#' 
#' @param support an \eqn{(n \times p)} matrix of fixed support where atoms are stacked as rows, or length-\eqn{n} vector.
#' @param measures a length-\eqn{K} list whose elements are \eqn{(n_i \times p)} matrix or length-\eqn{n_i} vector of atoms 
#' @param marginals marginal distribution for empirical measures; \code{NULL} for uniform weights or it should be a length-\eqn{K} list whose elements are length-\eqn{n_i} vectors that sum to 1.
#' @param weights weights for each individual measures; if \code{NULL}, it is set uniformly as \code{rep(1/K, K)} or a vector of length \eqn{K} should be provided.
#' @param p \eqn{p}-norm for computing distances.
#' @param maxiter maximum number of iterations as stopping criterion.
#' @param abstol threhold for stopping via incremental change of weights.
#' @param print.progress a logical; \code{TRUE} to print iteration.
#' @param useR if \code{TRUE}, we use \pkg{R} version of code in conjunction with partial \pkg{C++} while \code{FALSE} involves pure \pkg{C++} implementation.
#' 
#' @return a length-\eqn{n} vector of weights for atoms from \code{support}.
#' 
#' @examples 
#' ## simple example
#' #  empirical measures are created from N(0,I) in R^2
#' measures = list()
#' for (k in 1:10){
#'    ni = round(stats::runif(1, min=0.05, max=1)*100)
#'    measures[[k]] = matrix(rnorm(ni*2, mean=-1),ncol=2)
#' }
#' for (k in 11:20){
#'    ni = round(stats::runif(1, min=0.05, max=1)*100)
#'    measures[[k]] = matrix(rnorm(ni*2, mean=1),ncol=2)
#' }
#' 
#' # fixed support of 100 atoms 
#' # where domain variables are drawn from U(-3,3)
#' support = matrix(rnorm(100*2, sd=2), ncol=2)
#' 
#' # run
#' out1 = BaryProximal(support, measures, maxiter=500, useR=TRUE)
#' out2 = BaryProximal(support, measures, maxiter=500, useR=FALSE)
#'
#' par(mfrow=c(2,2))
#' plot(as.vector(out1))
#' plot(as.vector(out2))
#' 
#' p1 = as.vector(out1)
#' p2 = as.vector(out2)
#' id1 = sample(1:nrow(support), 50, prob=p1, replace=TRUE)
#' id2 = sample(1:nrow(support), 50, prob=p2, replace=TRUE)
#' plot(support[id1,])
#' plot(support[id2,])
#' 
#' 
#' @export
BaryProximal <- function(support, measures, marginals=NULL, weights=NULL, p=2, maxiter=496, abstol=1e-5, print.progress=TRUE, useR=FALSE){
  ##################################################
  # Preprocessing
  # 1. support
  if (!(is.vector(support)||is.matrix(support))){
    stop("* BaryProximal : 'support' should be either a vector or a matrix.")
  }
  # 2. measures
  if (!is.list(measures)){
    stop("* BaryProximal : 'measures' should be a list containing vectors/matrices.")
  }
  # 3. marginals : later
  # 4. weights   : later
  # 5. p
  myp = as.double(p)
  if ((length(myp)>1)||(myp<=0)){
    stop("* BaryProximal : 'p' should be a nonnegative real number.")
  }
  # 6. other ones
  myiter = round(maxiter)
  mytol  = max(double(abstol), sqrt(.Machine$double.eps))
  
  ##################################################
  # Preliminary Computation
  # 1. compute pairwise distances
  listdXY = extra_dist2list(support, measures, myp)
  # 2. compute cost matrices
  listcXY = list()
  for (i in 1:length(listdXY)){
    listcXY[[i]] = (listdXY[[i]]^myp)
  }
  # 3. check the marginals
  if ((length(marginals)==0)&&(is.null(marginals))){
    marginals = list()
    for (i in 1:length(listdXY)){
      ni = ncol(listcXY[[i]])
      marginals[[i]] = rep(1/ni, ni)
    }
  } else {
    if (length(marginals)!=length(listdXY)){
      stop("* BaryProximal : length of 'marginals' does not match with that of 'measures'.")
    }
    for (i in 1:length(listdXY)){
      if (length(marginals[[i]])!=ncol(listcXY[[i]])){
        stop(paste0("* BaryProximal : ",i,"-th marginal distribution does not have a matching length to the number of ",i,"-th atoms."))
      }
      if (sum(marginals[[i]])!=1){
        marginals[[i]] = marginals[[i]]/sum(marginals[[i]])
      }
    }
  }
  # 4. weights
  if ((length(weights)==0)&&(is.null(weights))){
    weights = rep(1/length(listdXY), length(listdXY))
  } 
  if (length(weights)!=length(listdXY)){
    stop(paste0("* BaryProximal : length of 'weights' is not matching to length of 'measures' ",length(listdXY)))
  }
  
  
  ##################################################
  # Main Computation
  # extra parameter
  N = nrow(listcXY[[1]])
  myiter = round(maxiter)
  mytol  = max(double(abstol), sqrt(.Machine$double.eps))
  # compute lambda
  lambda = 0.1
  for (k in 1:length(listcXY)){
    lambda = max(c(lambda, max(listcXY[[k]])/200))
  }
  # C++/R
  if (useR){
    # R version of iterations
    output = R_Barycenter_Proximal(N, myiter, mytol, listcXY, marginals, weights, lambda, as.logical(print.progress))
  } else {
    # C++ version of iterations
    output = C_Barycenter_Proximal(N, myiter, mytol, listcXY, marginals, weights, lambda, as.logical(print.progress))
  }
  
  ##################################################
  # Return
  return(as.vector(output))
}


#   -----------------------------------------------------------------------
#' @keywords internal
R_Barycenter_Proximal <- function(N, myiter, mytol, listcXY, marginals, weights, lambda, printer){
  # preprocess
  atil = rep(1/N, N)
  ahat = rep(1/N, N)
  t_0 = 2
  K   = length(listcXY)

  for (i in 1:myiter){
    a     = (1- 1/beta)*ahat + (1/beta)*atil
    beta  = (i+1)/2
    alpha = rep(0,N)
    for (k in 1:K){
      alpha <- alpha + weights[k]*extra_smooth_Subgradient(a, marginals[[k]], listcXY[[k]], lambda, 1000, 0.0001)  
    }
    atil = atil*exp(-(t_0)*beta*alpha)
    atil = atil/base::sum(atil)
    ahat = (1-1/beta)*ahat + (1/beta) * atil
    if (any(is.na(ahat))){
      return(a)
    }
    if (sqrt(sum((ahat-a)^2)) < mytol){
      return(ahat)
    }
    if (printer){
      print(paste0("* BaryProximal (R): iteration ",i,"/",myiter," complete.."))
    }
  }
  
  # run 
  return(ahat)
}
