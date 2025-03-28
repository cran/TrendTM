#' It performs the factorization for a fixed rank k and a temporal structure with a fixed tau
#'
#' @param Data_Series the data matrix with d rows and n columns containing the d temporal series with size n.
#' @param k the fixed rank of X. Default is 2.
#' @param tau the fixed value for tau . Default is \code{floor(n/2)}.
#' @param struct_temp a name indicating the temporal structure. Could be \code{none}, \code{periodic} or \code{smooth}. Default is \code{none}.
#' @param type_soft the option \code{type} of the function softImpute. Default is \code{als}.
#'
#' @return A list containing
#' \itemize{
#' \item \code{M_est} the estimation of M.
#' \item \code{U_est} the component U of the decomposition of \code{M_est}.
#' \item \code{V_est} the component V of the decomposition of \code{M_est}.
#' \item \code{contrast} the Frobenius norm of X-M_est.
#' }
#'
#' @export



FM_kt <- function(Data_Series, k = 2, tau = floor(n / 2), struct_temp = "none", type_soft = "als") {
  d <- nrow(Data_Series)
  n <- ncol(Data_Series)

  #####
  Lambda <- c()
  if (struct_temp == "none") {
    LambdaP <- c()
    Lambda <- diag(n)
    LambdaP <- Lambda
    X_tilde <- Data_Series %*% LambdaP
  }

  if (struct_temp == "periodic") {
    LambdaP <- c()
    p <- n / tau
    A <- matrix(1, ncol = p, nrow = 1)
    B <- diag(tau)
    Lambda <- kronecker(A, B)
    LambdaP <- kronecker(t(A), B) * (1 / p)
    X_tilde <- Data_Series %*% LambdaP
  }


  if (struct_temp == "smooth") {
    LambdaP <- c()
    times.eval <- seq(1 / n, 1, by = 1 / n)
    fbasis_obj <- fda::create.fourier.basis(rangeval = c(0, 1), nbasis = tau, period = 1)
    fbasis_evals <- fda::eval.basis(times.eval, fbasis_obj)
    Lambda <- t(fbasis_evals)
    LambdaP <- t(Lambda) / n
    X_tilde <- Data_Series %*% LambdaP
  }

  #####
  res_FM <- suppressWarnings(softImpute::softImpute(X_tilde, rank.max = k, lambda = 0, type_soft))

  if (k >= 2) {
    U_est <- (res_FM$u) %*% (diag(res_FM$d))
    V_est <- t(res_FM$v)
  }

  if (k == 1) {
    U_est <- as.matrix((res_FM$u) * res_FM$d)
    V_est <- t(as.matrix(res_FM$v))
  }

  M_tilde_est <- U_est %*% V_est
  M_est <- M_tilde_est %*% Lambda

  contrast <- c()
  contrast <- norm(Data_Series - M_est, type = "F")^2

  return(list(M_est = M_est, U_est = U_est, V_est = V_est, contrast = contrast))
}
