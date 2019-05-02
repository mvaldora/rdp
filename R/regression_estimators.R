############Regression estimators ###################################################
#' Computes a regression estimator of the mean of a random variable Y based on least squares
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept
#' @examples 
#'   vars <- rsamplemissing(n, alpha, beta)
#'   x <- vars$x
#'   y <- vars$y
#'   r(x, y)
#' @return The estimated mean of Y
#' @export

r <- function(x, y) {
  n <- length(y)
  a <- 1 - is.na(y)
  y[is.na(y)] <- 0
  yn <- y[a > 0]
  xn <- x[a > 0, ]
  ajuste <- lm.fit(cbind(1, xn), yn)
  betagorro <- ajuste$coefficients
  yt <- cbind(1, x) %*% betagorro
  mean(yt)
}
#' Computes a regression estimator of a quantile of a random variable Y with the robust method proposed in Sued and Yohai (2011)
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept
#' #' @examples 
#'   vars <- rsamplemissing(n, alpha, beta)
#'   x <- vars$x
#'   y <- vars$y
#'   a <- vars$a
#'   sy(x, y)$quantile
#' @return A list; element quantile is the estimated quantile of Y 
#' @export
sy <- function(x, y, qq=0.5)
{
  n <- length(y)
  a <- 1 - is.na(y)
  cara <- sum(a)
  yn <- y[!is.na(y)]
  xn <- x[!is.na(y), ]
  ycompleted <- y
  ycompleted[is.na(y)] <- 0
  ajustereg <-
    lmRob(yn ~ xn, control = lmRob.control(
      efficiency = 0.85,
      weight = c("Bisquare", "Bisquare")
    ))
  bet <- ajustereg$coefficients
  yt <- cbind(1, x) %*% bet
  res1 <- ycompleted - yt
  res <- res1[a > 0]
  v <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n)
  {
    v[i, ] <- (yt[i] + res1) * a
  }
  pesos <- rep(1 / (n * cara), n * cara)
  ans <- list()
  ans$pesos <- pesos
  pseudoyes <- c(v)
  pseudoyessinmissing <- v[, a == 1]
  ans$pseudoyes <- pseudoyes
  ans$pseudoyesmatrix <- v
  ans$pseudoyessinmissing <- pseudoyessinmissing
  #funcauxreg<-function(t,pseudosi) 1/cara*sum((pseudosi<=t)*a)
  funcreg <- function(t) {
    funcauxreg <- function(pseudosi)
      1 / cara * sum((pseudosi <= t) * a)
    apply(v, 1, funcauxreg)
  }
  ans$reg <- funcreg
  funcauxF <- function(t)
    mean(funcreg(t))
  ans$F <- Vectorize(funcauxF)
  ans$quantile <- discretequantile(qq, pseudoyessinmissing, pesos)
  return(ans)
}


#' Computes a regression estimator of a quantile of a random variable Y with the method proposed in Zhang et al (2011)
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept
#' @examples 
#'   vars <- rsamplemissing(n, alpha, beta)
#'   x <- vars$x
#'   y <- vars$y
#'   a <- vars$a
#'   rz(x, y)
#' @return The estimated quantile of Y 
#' @export

rz <- function(x, y, ps = NULL, qq = 0.5) {
  a <- 1 - is.na(y)
  cara <- sum(a)
  ycompleted <- y
  ycompleted[is.na(y)] <- 0
  if (is.null(ps))
    pred <- pigorro(x, a)
  else
    pred <- ps
  ajuste <- lm(y[a == 1] ~ x[a == 1, ])
  betagorro <- ajuste$coefficients
  sigmagorro <- sqrt(sum((y[a == 1] - ajuste$fitted) ^ 2 / (cara - 3)))
  vere <- cbind(1, x) %*% betagorro
  Fy <- function(t) {
    mean(pnorm((t - vere) / sigmagorro)) - qq
  }
  l1 <- min(y[a == 1]) - 1
  l2 <- max(y[a == 1])
  uniroot(Fy, c(l1, l2))$root
}
#' Computes a regression estimator of a quantile of a random variable Y with the robust method based on the method proposed Zhang et al (2011) described in Sued, Valdora and Yohai (2019)
#' @param y A vector of outcomes (may contain NAs)
#' @param x A matrix of covariates without an intercept
#' @examples 
#'   vars <- rsamplemissing(n, alpha, beta)
#'   x <- vars$x
#'   y <- vars$y
#'   a <- vars$a
#'   rrz(x, y)
#' @return The estimated quantile of Y 
#' @export
rrz <- function(x, y, ps = NULL, qq = 0.5) {
  a <- 1 - is.na(y)
  cara <- sum(a)
  if (is.null(ps))
    pred <- pigorro(x, a)
  else
    pred <- ps
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
  Fy <- function(t) {
    mean(pnorm((t - vere) / sigmagorro)) - qq
  }
  l1 <- min(y[a == 1]) - 1
  l2 <- max(y[a == 1])
  uniroot(Fy, c(l1, l2))$root
}





#example
if (FALSE) {
  vars <- rsamplemissing(n, alpha, beta)
  x <- vars$x
  y <- vars$y
  a <- vars$a
  reg(x, y)
  syflor(x, y)
  sy(x, y)$quantile
  rz(x, y)
  rrz(x, y)
}
