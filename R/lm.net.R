#' Title
#'
#' @param X 
#' @param Y 
#' @param Adj 
#' @param lambda 
#' @param gamma 
#'
#' @return A \code{list} object comprising:
#' \item{alpha} Vector of network cohesion.
#' \item{beta} Regression coefficients.
#' @export
#' 
#'
#' @examples
#' 
lm.net <- function(X, Y, Adj, lambda, gamma) {
  n <- nrow(X)
  p <- ncol(X)
  Y <- matrix(Y, ncol = 1)
  D <- diag(rowSums(Adj))
  L <- D - Adj

  Omega <- lambda * (L + gamma * diag(rep(1, n)))
  X.tilde <- cbind(diag(rep(1, n)), X)
  M <- matrix(0, n + p, n + p)
  M[1:n, 1:n] <- Omega
  alpha.beta <- solve(M + t(X.tilde) %*% X.tilde, t(X.tilde) %*% Y)
  alpha <- alpha.beta[1:n]
  beta <- alpha.beta[-(1:n)]
  return(list(alpha = alpha, beta = beta))
}
