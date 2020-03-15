#' Sinkhorn Distance's Subgradient
#' 
#' basic is for fixed point iteration
#' Bregman
#' 
#' @examples 
#' ## create two small datasets from bivariate normal
#' m = 20
#' n = 10
#' X = matrix(rnorm(m*2),ncol=2) # m obs. for X
#' Y = matrix(rnorm(n*2),ncol=2) # n obs. for Y
#' 
#' ## compute cross-distance between X and Y
#' dXY = array(0,c(m,n))
#' for (i in 1:m){
#'   vx = as.vector(X[i,])
#'   for (j in 1:n){
#'     vy  = as.vector(Y[j,])
#'     dXY[i,j] = sqrt(sum((vx-vy)^2))
#'   }
#' }
#' 
#' ## use several lambdas and compare with Wasserstein distance
#' outS1 = test_sub(dXY, lambda=0.01)
#' outS2 = test_sub(dXY, lambda=0.05)
#' outS3 = test_sub(dXY, lambda=0.1)
#' @export
test_sub <- function(dxy, p=1, wx=NULL, wy=NULL, lambda=0.1, maxiter=496, abstol=1e-6){
  ##################################################
  # Preprocessing
  # 1. dxy
  dxy = check_dxy(dxy, "Sinkhorn")
  m   = nrow(dxy)
  n   = ncol(dxy)
  # 2. p
  if ((length(p)>1)||(p<=0)){
    stop("* Sinkhorn : 'p' should be a nonnegative real number.")
  }
  # 3. wx and wy
  wx = check_weight(wx, m, "Sinkhorn")
  wy = check_weight(wy, n, "Sinkhorn")
  
  return(cpp_Sinkhorn_Subgradient(wx,wy,(dxy^p),1/lambda))
}