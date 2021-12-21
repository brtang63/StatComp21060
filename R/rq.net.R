#' @useDynLib StatComp21060
#' @importFrom Rcpp sourceCpp
NULL
#> NULL

#' @title Fitting Quantile Regression Models for Network-linked Data
#'
#' @param X Design matrix of predictors.
#' @param Y Response vector.
#' @param Adj Adjacency matrix.
#' @param lambda Regularization parameter.
#' @param tau The quantile to be estimated.
#'
#' @import CVXR
#'
#' @return A \code{list} object comprising:
#' \item{alpha}{Vector of network cohesion.} 
#' \item{beta}{Regression coefficients.} 
#'
#' @export
#'
#' @examples
#' V <- 100
#' K <- 3
#' p <- 3
#' data <- generate.data(V, K, p)
#' rq.net(data$X, data$Y, data$Adj)
rq.net <- function(X, Y, Adj, lambda = 0.1, tau = 0.5) {
  n <- nrow(X)
  p <- ncol(X)
  Y <- matrix(Y, ncol = 1)
  D <- diag(rowSums(Adj))
  L <- D - Adj

  # optimization using CVXR
  beta <- Variable(p)
  alpha <- Variable(n)
  quant_loss <- function(u, t) {
    0.5 * abs(u) + (t - 0.5) * u
  }
  obj <- sum(quant_loss(Y - alpha - X %*% beta, t = tau)) + lambda * quad_form(alpha, L)
  prob <- Problem(Minimize(obj))
  solution <- solve(prob)
  alpha <- solution$getValue(alpha)
  beta <- solution$getValue(beta)
  # print(solution$status)

  return(list(alpha = alpha, beta = beta))
}
