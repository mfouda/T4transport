#' Regularized Barycenter of Empirical Measures with Fixed Support via Sinkhorn Distance
#' 
#' 
#' @param support an \eqn{(n \times p)} matrix of fixed support where atoms are stacked as rows, or length-\eqn{n} vector.
#' @param measures a length-\eqn{K} list whose elements are \eqn{(n_i \times p)} matrix or length-\eqn{n_i} vector of atoms 
#' @param marginals marginal distribution for empirical measures; \code{NULL} for uniform weights or it should be a length-\eqn{K} list whose elements are length-\eqn{n_i} vectors that sum to 1.
#' @param weights weights for each individual measures; if \code{NULL}, it is set uniformly as \code{rep(1/K, K)} or a vector of length \eqn{K} should be provided.
#' @param lambda regularization parameter. If too small, it will lead to earlier break due to numerical precision issue.
#' @param p \eqn{p}-norm for computing distances.
#' @param maxiter maximum number of iterations as stopping criterion.
#' @param abstol threhold for stopping via incremental change of weights.
#' @param print.progress a logical; \code{TRUE} to print iteration.
#' @param useR if \code{TRUE}, we use \pkg{R} version of code in conjunction with partial \pkg{C++} while \code{FALSE} involves pure \pkg{C++} implementation.
#' 
#' @return a length-\eqn{n} vector of weights for atoms from \code{support}.
#' 
#' @examples 
#' ## example : two sets of empirical measures from Gaussians
#' #            whose means are (1,1) and (-1,-1).
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
#' # where domain variables are drawn from N(0,9I)
#' support = matrix(rnorm(100*2, sd=3), ncol=2)
#' 
#' # try different regularization levels
#' out1 = BarySinkhorn(support, measures, maxiter=100, lambda=0.1)
#' out2 = BarySinkhorn(support, measures, maxiter=100, lambda=0.5)
#' out3 = BarySinkhorn(support, measures, maxiter=100, lambda=1)
#'
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(1:100, out1, "b", pty=19, ylim=c(0,0.15), main="lambda=0.1")
#' plot(1:100, out2, "b", pty=19, ylim=c(0,0.15), main="lambda=0.5")
#' plot(1:100, out3, "b", pty=19, ylim=c(0,0.15), main="lambda=1")
#' par(opar)
#' 
#' @export
BarySinkhorn <- function(support, measures, marginals=NULL, weights=NULL, lambda=0.1, p=2, maxiter=496, abstol=1e-5, print.progress=TRUE, useR=TRUE){
  ##################################################
  # Preprocessing
  # 1. support
  if (!(is.vector(support)||is.matrix(support))){
    stop("* BarySinkhorn : 'support' should be either a vector or a matrix.")
  }
  # 2. measures
  if (!is.list(measures)){
    stop("* BarySinkhorn : 'measures' should be a list containing vectors/matrices.")
  }
  # 3. marginals : later
  # 4. weights   : later
  # 5. p
  myp = as.double(p)
  if ((length(myp)>1)||(myp<=0)){
    stop("* BarySinkhorn : 'p' should be a nonnegative real number.")
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
      stop("* BarySinkhorn : length of 'marginals' does not match with that of 'measures'.")
    }
    for (i in 1:length(listdXY)){
      if (length(marginals[[i]])!=ncol(listcXY[[i]])){
        stop(paste0("* BarySinkhorn : ",i,"-th marginal distribution does not have a matching length to the number of ",i,"-th atoms."))
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
    stop(paste0("* BarySinkhorn : length of 'weights' is not matching to length of 'measures' ",length(listdXY)))
  }
  
  
  ##################################################
  # Main Computation
  # extra parameter
  N = nrow(listcXY[[1]])
  myiter = round(maxiter)
  mytol  = max(double(abstol), sqrt(.Machine$double.eps))
  # compute lambda
  lambda = as.double(lambda)
  # C++/R
  if (useR){
    # R version of iterations
    output = R_Barycenter_Sinkhorn(N, myiter, mytol, listcXY, marginals, weights, lambda, as.logical(print.progress))
  } else {
    # C++ version of iterations
    output = C_Barycenter_Sinkhorn(N, myiter, mytol, listcXY, marginals, weights, lambda, as.logical(print.progress))
  }
  
  ##################################################
  # Return
  return(as.vector(output))
}

#   -----------------------------------------------------------------------
#' @keywords internal
R_Barycenter_Sinkhorn <- function(N, myiter, mytol, listcXY, marginals, weights, lambda, printer){
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
      print(paste0("* BarySinkhorn (R): iteration ",i,"/",myiter," complete.."))
    }
  }
  
  # run 
  return(ahat)
}
