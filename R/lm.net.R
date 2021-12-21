#' @title Fitting Linear Regression Models for Network-linked Data
#'
#' @param X Design matrix of predictors.``
#' @param Y Response vector.
#' @param Adj Adjacency matrix.
#' @param lambda Regularization parameter.
#' @param gamma Tuning parameter controlling the degree of jittering.
#'
#' @return A \code{list} object comprising:
#' \item{alpha}{Vector of network cohesion.} 
#' \item{beta}{Regression coefficients.} 
#' @export
#'
#'
#' @examples
#' V <- 100
#' p <- 3
#' K <- 3
#' data <- generate.data(V, K, p)
#' lm.net(data$X, data$Y, data$Adj)
lm.net <- function(X, Y, Adj, lambda = 0.1, gamma = 0.01) {
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
