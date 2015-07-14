#' Conduct optimal matching of cases to controls using network flow theory.
#' 
#' @param x.case a numeric \code{vector}.
#' @param x.control a numeric \code{vector}.
#' @param n a ratio of case:control.
#' @param method a function, a registry entry, or a mnemonic string referencing the proximity measure passed to \code{proxy::dist}.
#' @param ... further arguments passed to \code{proxy::dist}.
#' @references http://dx.doi.org/10.2307/2290079
#' @export
ccmatch <- function(x.case, x.control, n = 1, method = "Euclidean", ...) {
    if (any(!is.numeric(x.case)) || any(!is.numeric(x.control)))
        stop("Only numeric values are allowed.")
    if (nrow(x.case) == 0)
        stop("No case found.")
    if(nrow(x.control) == 0)
        stop("No control found.")
    if (n > nrow(x.control) / nrow(x.case))
        stop("Ratio of case:control is less than n.")    
    
    x.dist <- proxy::dist(x.case, x.control, method = method, ...)
    ret <- list()
    ret$ccmatch <- .ccmatch(x.dist, n)
    ret$score <- sum(ret$ccmatch$Distance)
    ret$Case <- ret$ccmatch$Case
    ret$Control <- sort(as.vector(t(ret$ccmatch[,2:(n+1)])))
    return (ret)
}