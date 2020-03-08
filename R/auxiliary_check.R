## Auxiliary Functions for 'T4transport' package
#   (1) check_dxy    : no negative, Inf, NA values as 'matrix'
#   (2) check_weight : simplex vector


# (1) check_dxy -----------------------------------------------------------
#' @keywords internal
check_dxy <- function(dxy, fnname){
  # 1. matrix
  if ((!is.matrix(dxy))||(length(dim(dxy))!=2)){
    stop(paste0("* ",fnname," : input 'dxy' should be a 2-dimensional matrix/array."))
  }
  # 2. values
  if (any(dxy < 0)){
    stop(paste0("* ",fnname," : input 'dxy' should contain no negative values."))
  }
  if (any(is.na(dxy))){
    stop(paste0("* ",fnname," : input 'dxy' should contain no NA values."))
  }
  if (any(is.infinite(dxy))){
    stop(paste0("* ",fnname," : input 'dxy' should contain no Inf values."))
  }
  return(as.matrix(dxy))
}

# (2) check_weight --------------------------------------------------------
#' @keywords internal
check_weight <- function(wvec, n, fnname){
  wname = deparse(substitute(wvec))
  if ((length(wvec)==0)&&(is.null(wvec))){
    return(rep(1/n, n))
  } else {
    wvec = wvec/sum(wvec)
    if (length(wvec)!=n){
      stop(paste0("* ",fnname," : input ",wname," is not a vector of length ",n,", not matching with 'dxy'."))
    }
    if (any(wvec<0)){
      stop(paste0("* ",fnname," : input ",wname," should contain no negative values."))
    }
    if (any(is.na(wvec))){
      stop(paste0("* ",fnname," : input ",wname," should contain no NA values."))
    }
    if (any(is.infinite(wvec))){
      stop(paste0("* ",fnname," : input ",wname," should contain no Inf values."))
    }
    return(wvec)
  }
}