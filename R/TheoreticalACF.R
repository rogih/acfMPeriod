#' Compute the theoretical ACOVF or ACF for a VAR(1) process
#'
#' Compute the theoretical autocovariance or autocorrelation function for a VAR(1) process.
#' @param sig a numeric white noise covariance matrix
#' @param phi a numeric coefficient matrix for the autoregressive order
#' @param lagmax maximum lag at which to calculate the acovf or acf
#' @param s1 an interger containing the entry of the series in the vector series. if s1 == s2 the acovf or acf is calculated.
#' @param s2 an interger containing the entry of the series in the vector series. if s1 != s2 the cross-acovf or cross-acf is calculated.
#' @param type character string giving the type to be computed. Allowed values are "correlation" (the default) or "covariance".
#' @return a numeric matrix
#' @author Higor Cotta. Adapted from Tsay's book.
#' @references Multivariate Time Series Analysis: with R and Financial Applications by Ruey S. Tsay, Wiley, 2014.
#' @export
#' @examples
#' phi <- matrix(c(0.7, 0.5, 0, 0.7), 2, 2)
#' cov.mat <- matrix(c(1, 0, 0, 1), 2, 2)
#' lagmax <- 9
#' TheoreticalACF(cov.mat, phi, lagmax, 1, 1, "correlation")
#' TheoreticalACF(cov.mat, phi, lagmax, 1, 2, "correlation")
#' TheoreticalACF(cov.mat, phi, lagmax, 2, 1, "correlation")
#' TheoreticalACF(cov.mat, phi, lagmax, 2, 2, "correlation")
TheoreticalACF <- function(sig, phi, lagmax, s1, s2, type = c("correlation", "covariance")) {
  type <- match.arg(type)
  nvar <- dim(phi)[1]
  Id <- diag(rep(1, nvar^2))
  kron <- kronecker(phi, phi)
  vecSig <- c(sig)
  vecSig <- matrix(vecSig, nvar^2, 1)
  Dif <- Id - kron
  DifInv <- solve(Dif)
  gamma0 <- DifInv %*% vecSig
  gamma0 <- matrix(gamma0, nvar, nvar)
  acf.ret <- matrix(NA, 1, lagmax)
  if (type == "covariance") {
    gammaBef <- gamma0
    acf.ret[, 1] <- gamma0[s1, s2]
    for (i in 2:lagmax) {
      gammaAtual <- phi %*% gammaBef
      gammaBef <- gammaAtual
      acf.ret[, i] <- gammaAtual[s1, s2]
    }
  }
  else {
    D <- diag(sqrt(diag(gamma0)))
    Di <- solve(D)
    rho0 <- Di %*% gamma0 %*% Di
    gammaBef <- gamma0
    acf.ret[, 1] <- rho0[s1, s2]
    for (i in 2:lagmax) {
      gammaAtual <- phi %*% gammaBef
      gammaBef <- gammaAtual
      rhoAtual <- Di %*% gammaAtual %*% Di
      acf.ret[, i] <- rhoAtual[s1, s2]
    }
  }
  return(acf.ret)
}
