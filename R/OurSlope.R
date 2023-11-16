
#' It performs the slope heuristic for the selection of a penalty constant
#'
#' @param grille the ordered grid of potential values for the penalty constant
#' @param contrast the Frobenius norm of X-M_est for all the value of the grid \code{grille}
#' @param penalty the penalty calculated for each value of the grid \code{grille}
#'
#' @return \code{Model_Selected} the selected parameter
#'
#' @export

OurSlope <- function(contrast, grille, penalty) {

  DataForCa <- data.frame(model=paste("K=",grille), pen=penalty, complexity=grille, contrast=contrast)
  Res_capushe <- suppressWarnings(capushe::DDSE(DataForCa)@ModelHat)
  rg_max <- which.max(Res_capushe$number_plateau)

  Model_Selected=grille[Res_capushe$model_hat[Res_capushe$point_breaking[rg_max]]]
  return(Model_Selected)

}


