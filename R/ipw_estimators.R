##########################################################
#' computes an ipw estimator of the quantiles of a random variable Y
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept (either \code{x} or \code{ps} muste be provided)
#' @param ps The estimated propensity score, i.e. the probability of y being observed conditional to \code{x} (either \code{x} or \code{ps} muste be provided)
#' @param type mean1 or mean2 corresponding to two ways of computing the ipw estimators, dividing by different constants, analogous to IPW1 and IPW2 in Lunceford and Davidian (2008) respectively 
#' @param p The desired quantile; default is 0.5, the median.
#' @return A list with elements pred (the propensity score), weights (the weight received by each observation) and median (the estimated median of Y)
#' @examples
#' n <- 100
#' alpha <- c(0, 0.1, -1.1)
#' beta <- c(0, -3, 10)
#' vars <- rsamplemissingnorm(n, alpha, beta)
#' x <- vars$x
#' y <- vars$y
#' a <- vars$a
#' quantileipw(y, x = x)$quantile
#' @export
quantileipw <- function(y,
                      x = NULL,
                      ps = NULL,
                      type = "mean2",p=0.5)
{
  ans <- list()
  n <- length(y)
  a <- 1 - is.na(y)
  y[is.na(y)] <- 0
  if (is.null(ps) &
      is.null(x)) {
    print("ERROR: Either ps or x must be provided")
  }
  if (is.null(ps)) {
    pred <- pigorro(x, a)
  } else
    pred <- ps
  if (sum(pred > 1 |
          pred < 0) > 0) {
    print("ERROR: ps values should be between 0 and 1")
  }
  if (type == "mean1") {
    pesos <- (a / pred) / n
    ans$pesos <- (a / pred) / n
  } else if (type == "mean2") {
    pesos <- (a / pred) / sum(a / pred)
  }
  ans$pred <- pred
  ans$weights <- pesos
  ans$quantile <- discretequantile(p, y, pesos)
  ans
}
#' computes an ipw estimator of the mean of a random variable Y
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept (either \code{x} or \code{ps} muste be provided)
#' @param ps The estimated propensity score, i.e. the probability of y being observed conditional to \code{x} (either \code{x} or \code{ps} muste be provided)
#' @param type mean1 or mean2 corresponding to two ways of computing the ipw estimators, dividing by different constants, analogous to IPW1 and IPW2 in Lunceford and Davidian (2008) respectively 
#' @return A list with elements pred (the propensity score), weights (the weight received by each observation) and mean (the estimated mean of Y)
#' @examples
#' n <- 100
#' alpha <- c(0, 0.1, -1.1)
#' beta <- c(0, -3, 10)
#' vars <- rsamplemissingnorm(n, alpha, beta)
#' x <- vars$x
#' y <- vars$y
#' a <- vars$a
#' meanipw(y, x = x)$mean
#' @export
meanipw <- function(y,
                    ps = NULL,
                    x = NULL,
                    type = "mean2")
{
  ans <- list()
  n <- length(y)
  a <- 1 - is.na(y)
  y[is.na(y)] <- 0
  if (is.null(ps) &
      is.null(x)) {
    print("ERROR: Either ps or x must be provided")
  }
  if (is.null(ps)) {
    pred <- pigorro(x, a)
  } else
    pred <- ps
  if (sum(pred > 1 |
          pred < 0) > 0) {
    print("ERROR: ps values should be between 0 and 1")
  }
  ans <- list()
  if (type == "mean1") {
    dis <- mean(y * (a / pred))
    pesos <- (a / pred) / n
  } else if (type == "mean2") {
    dis <- sum(y * (a / pred)) / sum(a / pred)
    pesos <- (a / pred) / sum(a / pred)
  }
  ans$mean <- dis
  ans$pred - pred
  ans$weights <- pesos
  ans
}

