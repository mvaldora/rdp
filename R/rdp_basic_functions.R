######################################################################
#' This function computes exp(x) / (1 + exp(x)), the inverse of the logistic function log(p / (1 - p))
#'
#' @param x A real number
#' @return A real number between 0 and 1
#' @export
expit <- function(x) {
  exp(x) / (1 + exp(x))
}


######################################################################
#' computes probability of a=1 conditional to x using the logistic model with parameter alpha
#'
#' @param x A matrix of covariates without the intercept
#' @param alpha A vector of parameters
#' @return a vector representing the probability of a=1 conditional to x using the logistic model with parameter alpha (with intercept automatically added)
#' @export
pax <- function(x, alpha)
{
  XI <- c(1, x)
  expit(sum(XI %*% alpha))
}
#' quantile function of a discrete random variable Y that takes on values y with probabilities py
#'
#' @param p A number in (0,1) representing the desired quantile (p=0.5 correspons to the median).
#' @param y A vector representing the rank of a discrete random variable.
#' @param py A vector of with elements in (0,1) that add up to 1, representing the probabilities that Y=y
#' @param ordered Logical; indicates whether or not the vector y is sorted in ascending order
#' @param  A list with elements x (the matrix of covariates), y (outcomes) and a (the missingness indicator)
#' @param  Logical; if TRUE, the vector y is assumed ordered, if FALSE it is not.
#' @return The \code{p}-quantile of Y
#' @export 
discretequantile <- function(p, y, py, ordered = FALSE) {
  if (ordered == TRUE)
    yy.orde <- y
  else
    yy.orde <- sort(y)
  pesos.orde <- py[order(y)]
  acumu <- cumsum(pesos.orde)
  indice <- which(acumu >= p)[1]# da la primer acumulada donde supera
  yy.orde[indice]
}
##########################################################
#' computes the propensity score
#' @param x A matrix of covariates without an intercept
#' @param a A vector of 0 and 1; if the i-th entry of a equals 1, it indicates the i-th outcome is missing.
#' @return A vector indicating the estimated probabilities that a equals 1 conditional to x.
#' @export
pigorro <- function(x, a)
{
  ajusteok <- glm(a ~ x, family = binomial(link = logit))
  pigok <- ajusteok$fitted.values
  pigok
}
##########################################################
#' computes the propensity score using the spps model introduced in Molina et al 2018
#' @param x A matrix of covariates without an intercept
#' @param a A vector of 0 and 1; if the i-th entry of a equals 1, it indicates the i-th outcome is missing.
#' @return A vector indicating the estimated probabilities that a equals 1 conditional to x.
#' @export
pigorrospps <- function(x, a) {
  ajuste <-
    estimall(
      x = x,
      a = a,
      phi = expit,
      dphi = dexpit,
      phiinv = expitinv,
      maxit = 20,
      estim.method = "ML",
      method = "IPW"
    )
  pred <- ajuste$fitted
  pred
}

