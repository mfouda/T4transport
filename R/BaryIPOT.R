#' Wasserstein Barycenter with IPOT Algorithm
#' 
#' 
#' 
#' @examples 
#' ## let's try a simple example
#' mydata = list()
#' for (i in 1:2){
#'   ni = round(stats::runif(1, min=10, max=20))
#'   mydata[[i]] = matrix(rnorm(ni*2),ncol=2)
#' }
#' hey = BaryIPOT(mydata)
#' 
#' 
#' @export
BaryIPOT <- function(list.data, list.weights=NULL, p=1, L=1, beta=0.1, lambda=rep(1/length(list.data), length(list.data))){
  ##################################################
  # Preprocessing
  # 1. data should be converted
  if (!is.list(list.data)){
    stop("* BaryIPOT : 'list.data' should be a list of data.")
  }
  K = length(list.data)
  mydata = c()
  mysize = c()
  for (k in 1:K){
    datak = list.data[[k]]
    if (is.vector(datak)){
      mydata = rbind(mydata, matrix(datak,ncol=1))
      mysize = c(mysize, length(datak))
    } else if (is.matrix(datak)){
      mydata = rbind(mydata, datak)
      mysize = c(mysize, nrow(datak))
    } else {
      stop(paste0("* BaryIPOT : ",k,"-th element of 'list.data' is not a vector nor a matrix."))
    }
  }
  
  # 2. weights
  allsize = sum(mysize)
  if ((length(list.weights)==0)&&(is.null(list.weights))){
    list.weights = list()
    for (k in 1:K){
      list.weights[[k]] = rep(1/mysize[k], mysize[k])
    }
  } 
  if (!is.list(list.weights)){
    stop("* BaryIPOT : 'list.weights' should be a vector.")
  }
  for (k in 1:K){
    list.weights[[k]] = list.weights[[k]]/sum(list.weights[[k]])
  }
  myweights = create_common_probability(list.weights)
  # 3. p
  if ((length(p)>1)||(p<=0)){
    stop("* BaryIPOT : 'p' should be a nonnegative real number.")
  }
  myp = as.double(p)
  dxy = extra_distmat(mydata)
  cxy = (dxy^myp)
  # 4. beta
  mybeta = as.double(beta)
  if (mybeta <= 0){
    stop("* BaryIPOT : 'beta' should be a nonnegative real number.")
  }
  if ((max(cxy)/mybeta) > 200){
    print("* BaryIPOT : 'beta' may be too small, inducing numerical instability.")
  }
  # 5. lambda 
  if (length(lambda)!=K){
    stop(paste0("* BaryIPOT : 'lambda' should be a vector of length ",K,"."))
  }
  # 6.
  myL = round(L)
  if (myL < 1){
    stop("* BaryIPOT : 'L' should be a nonnegative integer.")
  }
  ##################################################
  # Main Computation
  maxiter = 496
  abstol  = 1e-6
  Rweights = list()
  for (k in 1:K){
    Rweights[[k]] = as.vector(myweights[,k])
  }
  return(IPOTWB.R(C, Rweights, beta=mybeta))
  
    
  # return(cpp_IPOT_WB(cxy, myweights, mybeta, lambda, myL, maxiter, abstol))
}


#' R version of IPOT-WB
#' 
#' @examples 
#' X = matrix(rnorm(100*2),ncol=2)
#' C = as.matrix(stats::dist(X))^2
#' probs = list()
#' for (k in 1:20){
#'   tmprobs = stats::runif(100)
#'   probs[[k]] = tmprobs/sum(tmprobs)
#' }
#' 
#' @export
IPOTWB.R <- function(C, probs, beta=0.1){
  maxiter = 100
  abstol  = 1e-5
  K = length(probs)
  n = nrow(C)
  G = exp(-C/beta)
  
  Gammaold = array(1,c(n,n,K))
  Gammanew = array(0,c(n,n,K))
  list.H = list()
  list.a = list()
  list.b = list()
  for (k in 1:K){
    list.b[[k]] = rep(1/n,n)
  }
  q = rep(1/n,n)
  for (t in 1:maxiter){
    for (k in 1:K){
      list.H[[k]] = G*Gammaold[,,k]
    }
    for (k in 1:K){
      list.a[[k]] = q/as.vector(list.H[[k]]%*%list.b[[k]])
      list.b[[k]] = probs[[k]]/as.vector(t(list.H[[k]])%*%list.a[[k]])
    }
    q = rep(0,n)
    for (k in 1:K){
      q = q + (1/K)*log(list.a[[k]]*as.vector(list.H[[k]]%*%list.b[[k]]))
    }
    q = exp(q)
    for (k in 1:K){
      Gammanew[,,k] = diag(list.a[[k]])%*%list.H[[k]]%*%diag(list.b[[k]])
    }
    incGamma = diffGamma(Gammaold, Gammanew)
    print(paste0("iteration ",t," with incGamma=",incGamma))
    Gammaold = Gammanew
    if (incGamma < abstol){
      break
    }
  }
  return(q)
}
#' @export
IPOTWB.C <- function(C, probs, beta=0.1){
  maxiter = 100
  abstol  = 1e-5
  K = length(probs)
  n = nrow(C)
  G = exp(-C/beta)
  
  myweights = c()
  for (k in 1:K){
    myweights = cbind(myweights, probs[[k]])
  }
  mylambda = rep(1/K, K)
  return(as.vector(cpp_IPOT_WB(C, myweights, beta, mylambda, 1, maxiter, abstol)$q))
}
#' @keywords internal
diffGamma <- function(Gold, Gnew){
  output = 0
  K = dim(Gold)[3]
  for (k in 1:K){
    output = output + norm(Gold[,,k]-Gnew[,,k],"f")
  }
  return(output)
}

## personal example 