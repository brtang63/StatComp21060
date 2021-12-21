#' @title Generate simulated dataset for regression with network-linked data
#'
#' @param V number of vertices
#' @param K number of communities
#' @param inter.p in-community link probability
#' @param between.p out-community link probability
#'
#' @return A \code{list} object comprising:
#' \item{X}{Design matrix of predictors.} 
#' \item{Y}{Response vector.} 
#' \item{Adj}{Adjacency matrix.} 
#'
#' @import igraph
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rbinom
#'
#' @examples
#' V <- 100
#' K <- 3
#' dataset <- generate.data(V, K)
#' str(dataset)
#' @export
#'
generate.data <- function(V, K, inter.p = 0.8, between.p = 0.2) {
  p.x <- 1
  beta <- 1
  block.g <- BlockGen(V, K, inter.p, between.p)
  true.alpha <- AlphaGen.from.Block(block.g$Block, s = 0.2, scale = 3)
  dd <- Regression.Gen.norm(true.alpha$alpha, p = p.x, beta = beta)
  g <- igraph::graph_from_data_frame(block.g$Edge, directed = FALSE)
  Adj <- igraph::as_adjacency_matrix(g)
  Adj <- as.matrix(Adj)

  new.index <- as.numeric(igraph::V(g)$name)

  X <- matrix(dd$X[new.index, ], ncol = p.x)
  Y <- dd$Y[new.index]
  return(list(X = X, Y = Y, Adj = Adj))
}




BlockGen <- function(V, K, inter.p = 0.8, between.p = 0.2) {
  nodes <- rep(1:V)
  membership <- sample(K, size = V, replace = TRUE)
  node1 <- node2 <- NULL
  for (i in 1:(V - 1)) {
    for (j in (i + 1):V) {
      if (membership[i] == membership[j]) {
        connect <- stats::rbinom(1, 1, inter.p)
        if (connect == 1) {
          node1 <- c(node1, i)
          node2 <- c(node2, j)
        }
      } else {
        connect <- stats::rbinom(1, 1, between.p)
        if (connect == 1) {
          node1 <- c(node1, i)
          node2 <- c(node2, j)
        }
      }
    }
  }
  EdgeList <- data.frame(node1 = node1, node2 = node2)
  return(list(Edge = EdgeList, Block = membership))
}

AlphaGen.from.Block <- function(Block, s = 4, scale = 1) {
  K <- length(unique(Block))
  V <- length(Block)
  # scale <- s
  centroids <- seq(-scale, scale, length.out = K)
  alpha <- rep(0, V)
  for (v in 1:V) {
    alpha[v] <- stats::rnorm(1, mean = centroids[Block[v]], sd = s)
  }
  #    centroids <- centroids - mean(alpha)
  #    alpha <- alpha - mean(alpha)
  #    centroids <- centroids/(max(alpha)-min(alpha))*2*scale
  #    alpha <- alpha/(max(alpha)-min(alpha))*2*scale
  shift <- mean(alpha)
  alpha <- alpha - shift
  centroids <- centroids - shift
  return(list(alpha = alpha, centroids = centroids))
}

Regression.Gen.norm <- function(alpha, p, beta) {
  beta <- matrix(beta, ncol = 1)
  n <- length(alpha)
  X <- matrix(stats::runif(n * p, min = -2, max = 2), ncol = p)
  X <- scale(X, center = TRUE, scale = FALSE)
  EY <- X %*% beta + alpha
  Y <- EY + stats::rnorm(n)
  return(list(X = X, Y = Y, alpha = alpha, beta = beta, EY = EY))
}
