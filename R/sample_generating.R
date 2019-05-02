#' generates a sample with misssing values (x,a,y) that follows the distribution described in section  5 with normal errors
#'
#' @param n The sample size
#' @param alpha A vector of parameters corresponding to the logsitic model that relates \code{a} and \code{x}
#' @param beta A vector of parameters corresponding to the linear model that relates \code{y} and \code{x}
#' @return A list with elements x (the matrix of covariates), y (outcomes) and a (the missingness indicator)
#' @export
rsamplemissingnorm <- function(n, alpha, beta)
{
  p <- length(alpha) - 1
  x <- rmvnorm(n, mean = rep(0, p), sigma = diag(rep(1, p)))
  pcondax <- rep(0, n)
  for (i in 1:n)
  {
    pcondax[i] <- pax(x[i, ], alpha)
  }
  uni1 <- runif(n)
  a <- 1 * (uni1 < pcondax)
  u <- rnorm(n)
  XI <- cbind(1, x)
  y <- XI %*% beta + u
  y[a == 0] <- NA
  ans <- list()
  ans$x <- x
  ans$a <- a
  ans$y <- y
  ans
}
#' generates a sample with misssing values (x,a,y) that follows the distribution described in section  5 with errors following a t distribution
#'
#' @param n The sample size
#' @param alpha A vector of parameters corresponding to the logsitic model that relates \code{a} and \code{x}
#' @param beta A vector of parameters corresponding to the linear model that relates \code{y} and \code{x}
#' @param df The degrees of freedom of the t distribution of the errors
#' @return A list with elements x (the matrix of covariates), y (outcomes) and a (the missingness indicator)
#' @examples
#' n <- 100
#' alpha <- c(0, 0.1, -1.1)
#' beta <- c(0, -3, 10)
#' vars <- rsamplemissingnorm(n, alpha, beta)
#' @export
rsamplemissingt <- function(n, alpha, beta, df)
{
  p <- length(alpha) - 1
  x <- rmvnorm(n, mean = rep(0, p), sigma = diag(rep(1, p)))
  pcondax <- rep(0, n)
  for (i in 1:n)
  {
    pcondax[i] <- pax(x[i, ], alpha)
  }
  uni1 <- runif(n)
  a <- 1 * (uni1 < pcondax)
  u <- rt(n, df)
  XI <- cbind(1, x)
  y <- XI %*% beta + u
  y[a == 0] <- NA
  ans <- list()
  ans$x <- x
  ans$a <- a
  ans$y <- y
  ans
}
#' generates a sample with misssing values (x,a,y) that follows the distribution described in section  5 with errors following a cauchy distribution
#'
#' @param n The sample size
#' @param alpha A vector of parameters corresponding to the logsitic model that relates \code{a} and \code{x}
#' @param beta A vector of parameters corresponding to the linear model that relates \code{y} and \code{x}
#' @return A list with elements x (the matrix of covariates), y (outcomes) and a (the missingness indicator)
#' @export
rsamplemissingcauchy <- function(n, alpha, beta)
{
  p <- length(alpha) - 1
  x <- rmvnorm(n, mean = rep(0, p), sigma = diag(rep(1, p)))
  pcondax <- rep(0, n)
  for (i in 1:n)
  {
    pcondax[i] <- pax(x[i, ], alpha)
  }
  uni1 <- runif(n)
  a <- 1 * (uni1 < pcondax)
  u <- rcauchy(n)
  XI <- cbind(1, x)
  y <- XI %*% beta + u
  y[a == 0] <- NA
  ans <- list()
  ans$x <- x
  ans$a <- a
  ans$y <- y
  ans
}
