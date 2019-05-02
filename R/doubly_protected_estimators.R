#' computes a doubly protected estimator of the mean of a random variable Y with the robust method based on maximum likelihood methods
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept (either \code{x} or \code{ps} muste be provided)
#' @param ps The estimated propensity score, i.e. the probability of y being observed conditional to \code{x} (either \code{x} or \code{ps} muste be provided)
#' @examples 
#'   vars <- rsamplemissing(n, alpha, beta)
#'   x <- vars$x
#'   y <- vars$y
#'   a <- vars$a
#'   dp(x, y)
#' @return The estimated mean of Y 
#' @export

dp <- function(x, y, ps = NULL)
{
  a <- 1 - is.na(y)
  cara <- sum(a)
  ycompleted <- y
  ycompleted[is.na(y)] <- 0
  betagorro <- lm(y[a == 1] ~ x[a == 1, ])$coefficients
  vere <- cbind(1, x) %*% betagorro
  if (is.null(ps))
    pred <- pigorro(x, a)
  else
    pred <- ps
  mean(ycompleted * (a / pred)) - mean(vere * (a / pred)) + mean(vere)
}
#' Computes a robust doubly protected estimator of a quantile of a random variable Y with the robust method proposed in Sued, Valdora and Yohai (2019)
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept (either \code{x} or \code{ps} muste be provided)
#' @param ps The estimated propensity score, i.e. the probability of y being observed conditional to \code{x} (either \code{x} or \code{ps} muste be provided)
#' @param type 1 or 2; if type=2 the estimator is normalized as explained in Sued, Valdora and Yohai (2019)
#' @examples 
#'   vars <- rsamplemissing(n, alpha, beta)
#'   x <- vars$x
#'   y <- vars$y
#'   rdpn(x, y)
#' @return The estimated quantile of Y 
#' @export
rdpn <- function(x, y, ps = NULL, type = 2, qq = 0.5)
{
  a <- 1 - is.na(y)
  n <- length(y)
  cara <- sum(a)
  suedyohai <- sy(x, y , qq)
  vere <- suedyohai$pseudoyessinmissing
  if (is.null(ps))
    pred <- pigorro(x, a)
  else
    pred <- ps
  ans <- list()
  pseudoyes <- c(y[a == 1], vere)
  ans$pseudoyes <- pseudoyes
  pesospseudoyes <- matrix(0, nrow = n, ncol = n)
  if (type == 1) {
    cn <- n
  }
  if (type == 2) {
    cn <- sum(a / pred)
  }
  for (i in 1:n) {
    pesospseudoyes[, i] <- 1 / (n * cara) * a[i] - 1 / (cn * cara) * a / pred *
      a[i]
  }
  pesos <- c(a[a == 1] / pred[a == 1] / cn, pesospseudoyes[, a == 1])
  ans$weights <- pesos
  ans$quantile <- discretequantile(qq, pseudoyes, pesos)
  return(ans)
}


#' Computes a doubly protected estimator of a quantile of a random variable Y with the robust method proposed in Zhang et al. (2011)
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept (either \code{x} or \code{ps} muste be provided)
#' @param ps The estimated propensity score, i.e. the probability of y being observed conditional to \code{x} (either \code{x} or \code{ps} muste be provided)
#' @examples 
#'   vars <- rsamplemissing(n, alpha, beta)
#'   x <- vars$x
#'   y <- vars$y
#'   dpz(x, y)
#' @return The estimated quantile of Y 
#' @export
dpz <- function(x, y, ps = NULL, qq =0.5)
{
  a <- 1 - is.na(y)
  n <- length(y)
  cara <- sum(a)
  ycompleted <- y
  ycompleted[is.na(y)] <- 0
  ajuste <- lm(y[a == 1] ~ x[a == 1, ])
  betagorro <- ajuste$coefficients
  sigmagorro <- sqrt(sum((y[a == 1] - ajuste$fitted) ^ 2 / (cara - 3)))
  vere <- cbind(1, x) %*% betagorro
  if (is.null(ps))
    pred <- pigorro(x, a)
  else
    pred <- ps
  f <- function(t) {
    mean(pnorm((t - vere) / sigmagorro)) - qq + mean(a * ((ycompleted <= t) -
                                                             pnorm((t - vere) / sigmagorro)) / pred)
  }
  f <- Vectorize(f)
  l1 <- min(y[a == 1]) - 1
  l2 <- max(y[a == 1])
  uniroot(f, c(l1, l2))$root
}

#' Computes a robust doubly protected estimator of a quantile of a random variable Y with the method proposed in Zhang et al. (2011) robustified as explained in Sued, Valdora and Yohai (2019)
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept (either \code{x} or \code{ps} muste be provided)
#' @param ps The estimated propensity score, i.e. the probability of y being observed conditional to \code{x} (either \code{x} or \code{ps} muste be provided)
#' @examples 
#'   vars <- rsamplemissing(n, alpha, beta)
#'   x <- vars$x
#'   y <- vars$y
#'   pred <- pigorro(x, a)
#'   rdpz(x, y) #is the same as
#'   rdpz(x, y, ps = pred)
#' @return The estimated quantile of Y 
#' @export
rdpz <- function(x, y, ps = NULL, qq = 0.5)
{
  a <- 1 - is.na(y)
  n <- length(y)
  cara <- sum(a)
  ycompleted <- y
  ycompleted[is.na(y)] <- 0
  ajuste <-
    lmRob(y[a == 1] ~ x[a == 1, ], control = lmRob.control(
      efficiency = 0.85,
      weight = c("Bisquare", "Bisquare")
    ))
  betagorro <- ajuste$coefficients
  sigmagorro <- ajuste$scale
  vere <- cbind(1, x) %*% betagorro
  if (is.null(ps))
    pred <- pigorro(x, a)
  else
    pred <- ps
  f <- function(t) {
    mean(pnorm((t - vere) / sigmagorro)) - qq + mean(a * ((ycompleted <= t) -
                                                             pnorm((t - vere) / sigmagorro)) / pred)
  }
  f <- Vectorize(f)
  l1 <- min(y[a == 1]) - 1
  l2 <- max(y[a == 1])
  uniroot(f, c(l1, l2))$root
}

#examples
if (FALSE) {
  alpha <- c(0, 0.1, -1.1)
  beta <- c(0, -3, 10)
  n <- 100
  vars <- rsamplemissingnorm(n, alpha, beta)
  pred <- pigorro(x, a)
  x <- vars$x
  y <- vars$y
  a <- vars$a
  rdpn(x, y)$quantile
  dp(x, y) # is the same as
  dp(x, y, ps = pred)
  dpz(x, y)# is the same as
  dpz(x, y, ps = pred) # where pred is the propensity score (previously computed)
  rdpz(x, y)
  rdpz(x, y, ps = pred) # where pred is the propensity score (previously computed)
  rdpn(x, y)$quantile
  rdpn(x, y, ps = pred)$quantile # where pred is the propensity score (previously computed)
  dpz(x, y)  # is the same as
  dpz(x, y, ps = pred) # where pred is the propensity score (previously computed)
  rdpz(x, y) # is the same as
  rdpz(x, y, ps = pred) # where pred is the propensity score (previously computed)
}
