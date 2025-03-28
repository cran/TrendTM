#' Matrix Factorization for Multivariate Time Series Analysis
#'
#' It is the main function. It performs the factorization for a selected rank and a temporal structure with a selected tau if the selection is requested otherwise it is fixed
#'
#' @param Data_Series the data matrix with d rows and n columns containing the d temporal series with size n.
#' @param k_select a boolean indicating if the rank of the matrix Data_Series will be selected. Default is FALSE.
#' @param k_max the fixed rank of Data_Series if \code{k_select=FALSE}. The maximal value of the rank if \code{k_select=TRUE} (must be lower than the minimum between d and n). Default is 20.
#' @param struct_temp a name indicating the temporal structure. Could be \code{none}, \code{periodic} or \code{smooth}. Default is \code{none}.
#' @param tau_select a boolean indicating if the parameter tau will be selected. This can be possible only when \code{struct_temp=smooth}. Default is FALSE.
#' @param tau_max the fixed value for tau if \code{tau_select=FALSE}. The maximal value of tau if \code{tau_select=TRUE} (must be lower than n). Default is \code{floor(n/2)}.
#' @param type_soft the option \code{type} of the function softImpute. Default is \code{als}.
#'
#' @return A list containing
#' \itemize{
#' \item \code{k_est} the selected rank if \code{k_select==TRUE} or \code{k_max} if \code{k_select==FALSE}.
#' \item \code{tau_est} the selected tau if \code{tau_select==TRUE} or \code{tau_max} if \code{tau_select==FALSE}.
#' \item \code{U_est} the component U of the decomposition of the final estimator \code{M_est}.
#' \item \code{V_est} the component V of the decomposition of the final estimator \code{M_est}.
#' \item \code{M_est} the estimation of M.
#' \item \code{contrast} the Frobenius norm of Data_Series-M_est. This is a value when \code{k_select==FALSE} and \code{tau_select==FALSE}, a vector when \code{k_select==TRUE} or \code{tau_select==TRUE}, and a matrix when \code{k_select==TRUE} and \code{tau_select==TRUE} with \code{k_max} rows and \code{tau_max} columns.
#' }
#'
#' @details
#' The penalty constant(s) is(are) calibrated using the slope heuristic from package capushe. We adapt this heuristic as follows: the final dimension is the one correspind to the majority of the selected dimension for the considered different penalties.
#'
#' @examples
#' data(Data_Series)
#' result <- TrendTM(Data_Series, k_max = 3)
#' @export



TrendTM <- function(Data_Series, k_select = FALSE,
                    k_max = 20, struct_temp = "none", tau_select = FALSE,
                    tau_max = floor(n / 2), type_soft = "als") {

  d <- nrow(Data_Series)
  n <- ncol(Data_Series)
  s <- 4

  ######
  # Tests
  ## Data_Series
  if (is.matrix(Data_Series) == FALSE) {
    stop("Data_Series must be a matrix")
  }

  ## k
  if (is.numeric(k_max) == FALSE) {
    stop("k_max must be a number")
  }

  if (k_select != FALSE & k_select != TRUE) {
    stop("k_select must be a boolean TRUE or FALSE")
  }

  if (is.numeric(k_max) == TRUE & (k_max > min(d, n) | k_max < 1)) {
    stop("k_max must be a number between 1 and", min(d, n))
  }

  ## tau
  if (is.numeric(tau_max) == FALSE) {
    stop("tau_max must be a number")
  }

  if (is.numeric(tau_max) == TRUE & tau_max > n) {
    stop("tau_max must be a number inferior to", n)
  }

  if (tau_select != FALSE & tau_select != TRUE) {
    stop("tau_select must be a boolean TRUE or FALSE")
  }

  ## both
  if ((is.numeric(k_max) == TRUE) & (is.numeric(tau_max) == TRUE)) {
    if (k_max > tau_max) {
      stop("when no selection is required, k_max must be lower than tau_max")
    }
  }


  ## temporal structure
  if (struct_temp != "none" & struct_temp != "periodic" & struct_temp != "smooth") {
    stop("struct_temp must be none, periodic or smooth")
  }

  if (struct_temp == "smooth" & exists("tau_max") == FALSE) {
    stop("tau_max needs to be specified if there exists a temporal structure that is smooth")
  }

  if (struct_temp == "smooth" & is.null(tau_max) == TRUE) {
    stop("tau_max must be an integer if there exists a temporal structure that is smooth")
  }

  if (struct_temp != "smooth" & tau_select == TRUE) {
    stop("tau_select must be FALSE if there is no or periodic temporal structure")
  }

  if (tau_select == TRUE & struct_temp != "smooth") {
    stop("the selection of tau can only be done for a smooth temporal structure")
  }


  if (exists("tau_max")) {
    if (is.numeric(tau_max) == TRUE & struct_temp == "periodic") {
      if (n %% tau_max != 0) {
        stop("for periodic temporal structure, tau_max must be a number lower than", n, " and ", n, " must be a multiple of tau_max")
      }
    }
    if (is.numeric(tau_max) == TRUE & struct_temp == "smooth") {
      if (tau_max %% 2 == 0) {
        stop("for smooth temporal structure, tau_max must be an odd number lower than", n)
      }
    }
  }


  #######
  # No selection

  if (k_select == FALSE & tau_select == FALSE) {
    res <- FM_kt(Data_Series = Data_Series, k = k_max, tau = tau_max, struct_temp = struct_temp, type_soft = type_soft)
    k_est <- k_max
    tau_est <- n
    contrast <- res$contrast
  }

  ########
  # Fixed tau and selection of k
  if (k_select == TRUE & tau_select == FALSE) {
    grille_k <- 1:k_max
    contrast <- c()
    penalty <- c()
    penalty <- grille_k * (d + tau_max + s)

    contrast <- sapply(grille_k, function(j) FM_kt(Data_Series = Data_Series, k = j, tau = tau_max, struct_temp = struct_temp, type_soft = type_soft)$contrast)
    names(contrast) <- grille_k

    k_est <- OurSlope(contrast, grille_k, penalty)

    res <- FM_kt(Data_Series, k_est, tau_max, struct_temp, type_soft)
    tau_est <- tau_max
  }

  ########
  # Fixed k and selection of tau

  if (k_select == FALSE & tau_select == TRUE) {
    grille_tau <- c()
    if (k_max %% 2 == 0) {
      grille_tau <- seq(from = k_max + 1, to = tau_max, by = 2)
    } else {
      grille_tau <- seq(from = k_max + 2, to = tau_max, by = 2)
    }

    contrast <- c()
    penalty <- c()
    penalty <- k_max * (d + grille_tau + s)

    contrast <- sapply(grille_tau, function(j) FM_kt(Data_Series = Data_Series, k = k_max, tau = j, struct_temp = struct_temp, type_soft = type_soft)$contrast)
    names(contrast) <- grille_tau

    tau_est <- OurSlope(contrast, grille_tau, penalty)



    res <- FM_kt(Data_Series, k_max, tau_est, struct_temp, type_soft)
    k_est <- k_max
  }

  ########
  # Selection of k and tau


  if (k_select == TRUE & tau_select == TRUE) {
    grille_k <- 1:k_max
    grille_tau_min <- seq(from = 3, to = tau_max, by = 2)
    ntau_max <- length(grille_tau_min)
    contrast <- matrix(NA, ncol = ntau_max, nrow = k_max)
    contrast_k <- c()
    tau_est_k <- c()
    penalty_k <- c()

    for (i in 1:length(grille_k)) {
      k <- c()
      k <- grille_k[i]
      grille_tau <- c()
      if (k %% 2 == 0) {
        grille_tau <- seq(from = k + 1, to = tau_max, by = 2)
      } else {
        grille_tau <- seq(from = k + 2, to = tau_max, by = 2)
      }

      contrast.tau <- c()
      penalty.tau <- c()
      penalty.tau <- k * (d + grille_tau + s)

      contrast.tau <- sapply(grille_tau, function(j) FM_kt(Data_Series = Data_Series, k = k, tau = j, struct_temp = struct_temp, type_soft = type_soft)$contrast)
      contrast[i, grille_tau_min %in% grille_tau] <- contrast.tau


      tau_est_k[i] <- OurSlope(contrast.tau, grille_tau, penalty.tau)


      contrast_k[i] <- contrast.tau[grille_tau == tau_est_k[i]]
      penalty_k[i] <- penalty.tau[grille_tau == tau_est_k[i]]
    }


    rownames(contrast) <- grille_k
    colnames(contrast) <- grille_tau_min



    k_est <- OurSlope(contrast_k, grille_k, penalty_k)
    tau_est <- tau_est_k[which(k_est == grille_k)]
    res <- FM_kt(Data_Series, k_est, tau_est, struct_temp, type_soft)
  }

  return(list(k_est = k_est, tau_est = tau_est, M_est = res$M_est, U_est = res$U_est, V_est = res$V_est, contrast = contrast))
}
