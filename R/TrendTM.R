#' Matrix Factorization for Multivariate Time Series Analysis
#'
#'@param X the data matrix with d rows and n columns containing the d temporal series with size n.
#'
#' @param k.select a boolean indicating if the rank of the matrix X will be selected. Default is FALSE.
#' @param k.max the fixed rank of X if \code{k.select=FALSE}. The maximal value of the rank if \code{k.select=TRUE} (must be lower than the minimum between d and n). Default is 20.
#' @param struct.temp a name indicating the temporal structure. Could be \code{none}, \code{periodic} or \code{smooth}. Default is \code{none}.
#' @param tau.select a boolean indicating if the parameter tau will be selected. This can be possible only when \code{struct.temp=smooth}. Default is FALSE.
#' @param tau.max the fixed value for tau if \code{tau.select=FALSE}. The maximal value of tau if \code{tau.select=TRUE} (must be lower than n). Default is \code{floor(n/2)}.
#' @param type.soft the option \code{type} of the function softImpute. Default is \code{als}.
#'
#' @return A list containing
#' \itemize{
#' \item \code{k.est} the selected rank if \code{k.select==TRUE} or \code{k.max} if \code{k.select==FALSE}.
#' \item \code{tau.est} the selected tau if \code{tau.select==TRUE} or \code{tau.max} if \code{tau.select==FALSE}.
#' \item \code{U.est} the component U of the decomposition of the final estimator \code{M.est}.
#' \item \code{V.est} the component V of the decomposition of the final estimator \code{M.est}.
#' \item \code{M.est} the estimation of M.
#' \item \code{contrast} the Frobenius norm of X-M.est. This is a value when \code{k.select==FALSE} and \code{tau.select==FALSE}, a vector when \code{k.select==TRUE} or \code{tau.select==TRUE}, and a matrix when \code{k.select==TRUE} and \code{tau.select==TRUE} with \code{k.max} rows and \code{tau.max} columns.
#' }
#'
#' @details
#' The penalty constant(s) is(are) calibrated using the slope heuristic from package capushe. We adapt this heuristic as follows: the final dimension is the one correspind to the majority of the selected dimension for the considered different penalties.
#'
#' @examples
#' data(DataX)
#' k.max=3
#' result=TrendTM(X, k.max=k.max)
#' @export


#########################
### FONCTION TREND.TM ###
#########################

TrendTM <- function(X,k.select=FALSE, k.max= 20, struct.temp="none",tau.select=FALSE, tau.max=floor(n/2),type.soft="als"){

  d = nrow(X)
  n = ncol(X)
  s=4

  ######
  # Tests
    ## X
  if (is.matrix(X) == FALSE){
    stop('X must be a matrix')
  }

    ## k
  if (is.numeric(k.max) == FALSE){
    stop('k.max must be a number')
  }

  if (k.select != FALSE & k.select != TRUE){
    stop('k.select must be a boolean TRUE or FALSE')
  }

  if (is.numeric(k.max) == TRUE   & (k.max > min(d,n) | k.max <1)){
    stop('k.max must be a number between 1 and',min(d,n))
  }

    ## tau
  if (is.numeric(tau.max) == FALSE){
    stop('tau.max must be a number')
  }

  if (is.numeric(tau.max) == TRUE   & tau.max > n){
    stop('tau.max must be a number inferior to',n)
  }

  if (tau.select != FALSE & tau.select != TRUE){
    stop('tau.select must be a boolean TRUE or FALSE')
  }

    ## both
  if ((is.numeric(k.max) == TRUE) & (is.numeric(tau.max) == TRUE)){
      if (k.max > tau.max)
        stop('when no selection is required, k.max must be lower than tau.max')
  }


    ## temporal structure
  if (struct.temp != "none" & struct.temp != "periodic" & struct.temp != "smooth"){
    stop('struct.temp must be none, periodic or smooth')
  }

  if (struct.temp == "smooth" & exists("tau.max")==FALSE){
    stop('tau.max needs to be specified if there exists a temporal structure that is smooth')
  }

  if (struct.temp == "smooth" & is.null(tau.max)==TRUE ) {
    stop('tau.max must be an integer if there exists a temporal structure that is smooth')
  }

  if (struct.temp != "smooth" & tau.select==TRUE ) {
    stop('tau.select must be FALSE if there is no or periodic temporal structure')
  }

  if (tau.select == TRUE & struct.temp != "smooth"){
    stop('the selection of tau can only be done for a smooth temporal structure')
  }


  if (exists("tau.max")){
    if (is.numeric(tau.max) == TRUE & struct.temp=="periodic"){
      if (n%%tau.max !=0){
        stop('for periodic temporal structure, tau.max must be a number lower than',n,' and ',n,' must be a multiple of tau.max')
      }
    }
    if (is.numeric(tau.max) == TRUE & struct.temp=="smooth"){
      if ( tau.max%%2==0){
        stop('for smooth temporal structure, tau.max must be an odd number lower than',n)
      }
    }
  }


  #######
  # No selection

  if (k.select == FALSE & tau.select == FALSE){
    res <-  FM.kt(X=X, k=k.max,tau=tau.max,struct.temp=struct.temp,type.soft=type.soft)
    k.est <-  k.max
    tau.est <-  n
    contrast <- res$contrast
  }

  ########
  # Fixed tau and selection of k
  if (k.select == TRUE & tau.select == FALSE){

    grille_k <-  1:k.max
    contrast <- c()
    penalty <- c()
    penalty <- grille_k*(d+tau.max+s)

    contrast <- sapply(grille_k,function(j) FM.kt(X=X,k=j, tau= tau.max,struct.temp=struct.temp,type.soft=type.soft)$contrast)
    names(contrast) <- grille_k

    k.est=OurSlope(contrast, grille_k, penalty)

    res <-  FM.kt(X, k.est,tau.max,struct.temp, type.soft)
    tau.est <-  tau.max
  }

  ########
  # Fixed k and selection of tau

  if (k.select == FALSE & tau.select == TRUE){

    grille_tau <-  c()
    if (k.max%%2==0){
      grille_tau <-  seq(from = k.max+1, to = tau.max, by = 2)
      } else {grille_tau <-  seq(from = k.max+2, to = tau.max, by = 2)}

    contrast <- c()
    penalty <- c()
    penalty <- k.max*(d+grille_tau+s)

  contrast <- sapply(grille_tau,function(j) FM.kt(X=X,k=k.max, tau= j,struct.temp=struct.temp,type.soft=type.soft)$contrast)
  names(contrast) <- grille_tau

  tau.est=OurSlope(contrast, grille_tau, penalty)



    res <-  FM.kt(X, k.max,tau.est,struct.temp, type.soft)
    k.est <-  k.max
  }

  ########
  # Selection of k and tau


  if (k.select == TRUE & tau.select == TRUE){

    grille_k <-  1:k.max
    grille_tau_min <- seq(from = 3, to = tau.max, by = 2)
    ntau.max=length(grille_tau_min)
    contrast <- matrix(NA,ncol= ntau.max,nrow=k.max)
    contrast_k <- c()
    tau.est_k <- c()
    penalty_k <- c()

    for (i in 1:length(grille_k)){

      k <- c()
      k <- grille_k[i]
      grille_tau = c()
      if (k%%2==0){
        grille_tau <-  seq(from = k+1, to = tau.max, by = 2)
      } else {grille_tau <-  seq(from = k+2, to = tau.max, by = 2)
        }

      contrast.tau <- c()
      penalty.tau <- c()
      penalty.tau <- k*(d+grille_tau+s)

      contrast.tau <- sapply(grille_tau,function(j) FM.kt(X=X,k=k, tau= j,struct.temp=struct.temp,type.soft=type.soft)$contrast)
      contrast[i,grille_tau_min %in% grille_tau] <- contrast.tau


      tau.est_k[i]=OurSlope(contrast.tau, grille_tau, penalty.tau)


      contrast_k[i] <- contrast.tau[grille_tau==tau.est_k[i]]
      penalty_k[i] <- penalty.tau[grille_tau==tau.est_k[i]]
}


    rownames(contrast) <- grille_k
    colnames(contrast) <- grille_tau_min



    k.est=OurSlope(contrast_k, grille_k, penalty_k)
    tau.est = tau.est_k[which(k.est==grille_k)]
    res <-  FM.kt(X, k.est,tau.est,struct.temp, type.soft)
    }

  return(list(k.est=k.est,tau.est=tau.est,M.est = res$M.est, U.est = res$U.est, V.est=res$V.est, contrast=contrast))

}

# Used functionds

# Estimation of M for a fixed k and a fixed tau
FM.kt <- function(X,  k=2, tau=n, struct.temp="none",type.soft="als"){

  ####
  # Dimension of X
  d = nrow(X)
  n = ncol(X)

  #####
  Lambda=c()
  if (struct.temp=="none"){
    LambdaP=c()
    Lambda <-  diag(n)
    LambdaP <- Lambda
    X_tilde <-  X%*%LambdaP
  }

  if (struct.temp=="periodic"){
    LambdaP=c()
    p <- n/tau
    A <- matrix(1, ncol=p,nrow=1)
    B <- diag(tau)
    Lambda <- kronecker(A,B)
    LambdaP <- kronecker(t(A),B)*(1/p)
    X_tilde <-  X%*%LambdaP
  }


  if (struct.temp=="smooth"){
    LambdaP=c()
    times.eval <-seq(1/n,1,by=1/n)
    fbasis_obj= fda::create.fourier.basis(rangeval=c(0,1),nbasis=tau,period = 1)
    fbasis_evals <- fda::eval.basis(times.eval, fbasis_obj)
    Lambda <- t(fbasis_evals)
    LambdaP <- t(Lambda)/n
    X_tilde <-  X%*%LambdaP
  }

  #####
  res.FM = suppressWarnings(softImpute::softImpute(X_tilde,rank.max=k,lambda=0,type.soft))
  if (k>=2){
    U.est = (res.FM$u)%*%(diag(res.FM$d))
    V.est=t(res.FM$v)
  }

  if (k==1){
    U.est = as.matrix((res.FM$u)*res.FM$d)
    V.est = t(as.matrix(res.FM$v))
  }

  M_tilde_est = U.est%*%V.est
  M_est = M_tilde_est%*%Lambda

  contrast <- c()
  contrast <- norm(X-M_est,type="F")^2

  return(list(M.est=M_est, U.est=U.est, V.est=V.est, contrast=contrast))

}


# Adapted Slope heuristic
OurSlope <- function(contrast, grille, penalty){



  DataForCa=data.frame(model=paste("K=",grille),pen=penalty,complexity=grille,contrast=contrast)
  Res.capushe=suppressWarnings(capushe::DDSE(DataForCa)@ModelHat)
  rg.max=which.max(Res.capushe$number_plateau)

  Model.Selected=grille[Res.capushe$model_hat[Res.capushe$point_breaking[rg.max]]]
  return(Model.Selected)

}


