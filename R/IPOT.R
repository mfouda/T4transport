#' IPOT
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
#' ## use several betas and compare with Wasserstein distance
#' outLP = Wasserstein(dXY, method="lp")
#' outP1 = IPOT(dXY, beta=0.01)
#' outP2 = IPOT(dXY, beta=0.05)
#' outP3 = IPOT(dXY, beta=0.1)
#' 
#' ## visualize
#' #  show the transport plans and computed distance
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' image(outLP$plan,main=paste0("LP:",round(outLP$distance,3)))
#' image(outP1$plan,main=paste0("beta=0.01:",round(outP1$distance,3)))
#' image(outLP$plan,main=paste0("beta=0.05:",round(outP2$distance,3)))
#' image(outLP$plan,main=paste0("beta=0.10:",round(outP3$distance,3)))
#' par(opar)
#' 
#' @export
IPOT <- function(dxy, p=1, wx=NULL, wy=NULL, beta=0.1, L=1, maxiter=496, abstol=1e-6, print.progress=FALSE){
  ##################################################
  # Preprocessing
  # 1. dxy
  dxy = check_dxy(dxy, "IPOT")
  m   = nrow(dxy)
  n   = ncol(dxy)
  # 2. p
  if ((length(p)>1)||(p<=0)){
    stop("* IPOT : 'p' should be a nonnegative real number.")
  }
  # 3. wx and wy
  wx = check_weight(wx, m, "IPOT")
  wy = check_weight(wy, n, "IPOT")
  # 4. beta
  mybeta = as.double(beta)
  if (mybeta <= 0){
    stop("* IPOT : 'beta' should be a nonnegative real number.")
  }
  if ((max((dxy)^p)/beta) > 200){
    print("* IPOT : 'beta' may be too small, inducing numerical instability.")
  }
  # 6. other variables
  myiter = round(maxiter)
  mytol  = as.double(abstol)
  myL    = round(L)
  
  ##################################################
  # Main Computation
  if (is.infinite(p)){
    output = list()
    output$distance = max(dxy)
  } else {
    cxy = (dxy^p)
    output = cpp_IPOT(wx, wy, cxy, mybeta, myL, myiter, mytol, as.logical(print.progress))
    output$distance = (output$distance^(1/p))
  }
  return(output)
}