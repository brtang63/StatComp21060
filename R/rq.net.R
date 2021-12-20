#' Title
#'
#' @param X 
#' @param Y 
#' @param Adj 
#' @param lambda 
#' @param tau 
#' 
#' @import CVXR
#'
#' @return A \code{list} object comprising:
#' \item{alpha} Vector of network cohesion.
#' \item{beta} Regression coefficients.
#'
#' @export
#'
#' @examples
#' 
#' 
rq.net <- function(X, Y, Adj, lambda, tau = 0.5) {
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
