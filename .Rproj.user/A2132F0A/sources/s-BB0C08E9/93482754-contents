#' Wasserstein Distance between Empirical Measures
#' 
#' @param dxy an \eqn{(M\times N)} matrix of distances \eqn{d(x_m,y_n)}.
#' @param p an exponent for Wasserstein distance \eqn{\mathcal{W}_p}.
#' @param wx a length-\eqn{M} weight vector that sums to 1; if \code{NULL}, it automatically takes uniform weights \code{rep(1/M,M)}.
#' @param wy a length-\eqn{N} weight vector that sums to 1; if \code{NULL}, it automatically takes uniform weights \code{rep(1/N,N)}.
#' @param method name of an algorithm
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
#' ## compute algorithms for W_1 distances
#' out1 = Wasserstein(dXY, method="lp")
#' out2 = Wasserstein(dXY, method="networkflow")
#' out3 = Wasserstein(dXY, method="primaldual")
#' 
#' ## visualize
#' #  show the transport plan and computed distance
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' image(out1$plan,main=paste0("LP:",round(out1$distance,4)))
#' image(out2$plan,main=paste0("Network Flow:",round(out2$distance,4)))
#' image(out3$plan,main=paste0("Primal-Dual:",round(out3$distance,4)))
#' par(opar)
#' 
#' @export
Wasserstein <- function(dxy, p=1, wx=NULL, wy=NULL,
                        method=c("lp","networkflow","shortsimplex","revsimplex","primaldual")){
  ##################################################
  # Preprocessing
  # 1. dxy
  dxy = check_dxy(dxy, "Wasserstein")
  m   = nrow(dxy)
  n   = ncol(dxy)
  # 2. p
  if ((length(p)>1)||(p<=0)){
    stop("* Wasserstein : 'p' should be a nonnegative real number.")
  }
  # 3. wx and wy
  wx = check_weight(wx, m, "Wasserstein")
  wy = check_weight(wy, n, "Wasserstein")
  # 4. method
  mymethod = match.arg(method)

  ##################################################
  # Main Computation
  if (is.infinite(p)){
    output = list()
    output$distance = max(dxy)
  } else {
    cxy = (dxy^p)
    output = switch(mymethod,
                    "lp" = Wasserstein_lp(cxy, wx, wy),
                    "networkflow"  = Wasserstein_transport(cxy, wx, wy, mymethod),
                    "shortsimplex" = Wasserstein_transport(cxy, wx, wy, mymethod),
                    "revsimplex"   = Wasserstein_transport(cxy, wx, wy, mymethod),
                    "primaldual"   = Wasserstein_transport(cxy, wx, wy, mymethod))
    output$distance = (output$distance^(1/p))
  }

  ##################################################
  # Return
  return(output)
}


# auxiliary functions -----------------------------------------------------
#' @keywords internal
Wasserstein_lp <- function(cxy, wx, wy){
  m = nrow(cxy)
  n = ncol(cxy)
  
  c = as.vector(cxy)
  A = rbind(base::kronecker(matrix(1,nrow=1,ncol=n), diag(m)),
            base::kronecker(diag(n), matrix(1,nrow=1,ncol=m)))
  
  f.obj = c
  f.con = A
  f.dir = rep("==",n)
  f.rhs = c(rep(1/m,m),rep(1/n,n))
  f.sol = lp("min", f.obj, f.con, f.dir, f.rhs)
  
  gamma = matrix(f.sol$solution, nrow=m)
  value = sum(gamma*cxy)
  
  return(list(distance=value, plan=gamma))
}
#' @keywords internal
Wasserstein_transport <- function(cxy, wx, wy, mymethod){
  nx = nrow(cxy)
  ny = ncol(cxy)
  output = array(0,c(nx,ny))
  
  myplan = transport::transport(a=wx, b=wy, costm=cxy, method=mymethod)
  nplans = nrow(myplan)
  for (i in 1:nplans){
    id.from = myplan[i,1]
    id.to   = myplan[i,2]
    
    output[id.from,id.to] = myplan[i,3]
  }
  
  result = list()
  result$distance = sum(output*cxy)
  result$plan     = output
  return(result)
}
#' @keywords internal
Wasserstein_sinkhorn <- function(cxy, wx, wy, eps){
  maxiter = 100
  abstol  = 1e-6
  Gxy = exp(-cxy/eps)
  
  plan = cpp_sinkhorn(Gxy, wx, wy, maxiter, abstol)
  dist = sum(plan*cxy)
  return(list(distance=dist, plan=plan))
}