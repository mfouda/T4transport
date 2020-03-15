#' Proximal Point Algorithm 
#' 
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
#' ## use several betas for exact and inexact algorithms
#' outE1 = Proximal(dXY, beta=0.01, method="exact")
#' outE2 = Proximal(dXY, beta=0.05, method="exact")
#' outE3 = Proximal(dXY, beta=0.10, method="exact")
#' outI1 = Proximal(dXY, beta=0.01, method="inexact")
#' outI2 = Proximal(dXY, beta=0.05, method="inexact")
#' outI3 = Proximal(dXY, beta=0.10, method="inexact")
#' 
#' ## visualize
#' #  distances
#' dE1 = round(outE1$distance, 3) # exact distances
#' dE2 = round(outE2$distance, 3)
#' dE3 = round(outE3$distance, 3)
#' dI1 = round(outI1$distance, 3) # inexact distances
#' dI2 = round(outI2$distance, 3)
#' dI3 = round(outI3$distance, 3)
#' 
#' #  show the transport plans and computed distance
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3))
#' image(outE1$plan,main=paste0("beta=0.01:",dE1), xlab="Exact")
#' image(outE2$plan,main=paste0("beta=0.05:",dE2), xlab="Exact")
#' image(outE3$plan,main=paste0("beta=0.10:",dE3), xlab="Exact")
#' image(outI1$plan,main=paste0("beta=0.01:",dI1), xlab="Inexact")
#' image(outI2$plan,main=paste0("beta=0.05:",dI2), xlab="Inexact")
#' image(outI3$plan,main=paste0("beta=0.10:",dI3), xlab="Inexact")
#' par(opar)
#' 
#' @export
Proximal <- function(dxy, p=1, wx=NULL, wy=NULL, beta=0.1, L=1, 
                     method=c("inexact","exact"), maxiter=496, abstol=1e-6, print.progress=FALSE){
  ##################################################
  # Preprocessing
  # 1. dxy
  dxy = check_dxy(dxy, "Proximal")
  m   = nrow(dxy)
  n   = ncol(dxy)
  # 2. p
  if ((length(p)>1)||(p<=0)){
    stop("* Proximal : 'p' should be a nonnegative real number.")
  }
  # 3. wx and wy
  wx = check_weight(wx, m, "Proximal")
  wy = check_weight(wy, n, "Proximal")
  # 4. beta
  mybeta = as.double(beta)
  if (mybeta <= 0){
    stop("* Proximal : 'beta' should be a nonnegative real number.")
  }
  if ((max((dxy)^p)/beta) > 200){
    print("* Proximal : 'beta' may be too small, inducing numerical instability.")
  }
  # 6. other variables
  myiter = round(maxiter)
  mytol  = as.double(abstol)
  myL    = round(L)
  # 7. method=c("exact","inexact")
  if (missing(method)){
    mymethod = "inexact"
  } else {
    allmethod=c("exact","inexact")
    mymethod = match.arg(tolower(method), allmethod)  
  }
  
  ##################################################
  # Main Computation
  if (is.infinite(p)){
    output = list()
    output$distance = max(dxy)
  } else {
    cxy = (dxy^p)
    if (all(mymethod=="inexact")){
      output = cpp_IPOT(wx, wy, cxy, mybeta, myL, myiter, mytol, as.logical(print.progress))
    } else if (all(mymethod=="exact")){
      output = cpp_EPOT(wx, wy, cxy, mybeta, myL, myiter, mytol, as.logical(print.progress))
    }
    output$distance = (output$distance^(1/p))  
  }
  return(output)
}