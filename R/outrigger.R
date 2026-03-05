outrigger <- function(formula, data, xtest, degree=0, bandwidth=NULL, kernel="epan", lambda=NULL, folds=5, scoregranality=1, score_df=10){
  #' Performs outrigger regression at a point x
  #'
  #' @param formula formula (as in e.g. \code{lm}).
  #' @param data data frame.
  #' @param xtest datapoint(s) for local estimation.
  #' @param degree degree of polynomial.
  #' @param bandwidth bandwidth parameter. If NULL a rule-of-thumb bandwidth is selected by least-squares cross-validation for the standard local polynomial estimator.
  #' @param kernel kernel for local polyomial. Default is epanechnikov kernel.
  #' @param lambda orthogonalisation parameter. By default takes \eqn{\lambda = 4(1 + 0.5 \log(n))}.
  #' @param folds number of folds for cross-fitting. Default is 5.
  #' @param score_df number of degrees of freedom for score matching splines used for conditional score estimation. Input can be a numeric value (e.g. 10 on larger datasets, 6 on smaller datasets) or "cross-validation" or "cross-validation-quick", in which case score-matching-CV will be performed (may be computationally costly).
  #' @param scoregranality the number of bins to split covariate-space \eqn{|X-x|\leq\lambda h} for conditional score estimation.
  #'
  #' @return If x is a single point, return a list containing:
  #'   \describe{
  #'     \item{\code{prediction}}{The outrigger local polynomial estimator \eqn{\hat{f}(x)} at \eqn{x}.}
  #'     \item{\code{fitted_vector}}{For the local linear outrigger estimator, the full fitted vector at \eqn{x}}
  #'     \item{\code{standardlocpol_fitted_vector}}{The fitted vector for the standard local polynomial}
  #'     \item{\code{score_plot_metadata}}{Data used for \code{score_plotting} to plot the fitted conditional score functions}
  #'   }
  #'   If x is a vector of points, return a dataset with points x and associated outrigger fitted values.
  #' @importFrom np npregbw
  #' @importFrom stats quantile
  #' @export
  #'

  if (length(xtest)==1) {
    return(outrigger_onepoint(formula, data, xtest, degree, bandwidth, kernel, lambda, folds, scoregranality, score_df))
  } else if (length(xtest)>1) {
    ntest <- length(xtest)
    if (ntest*nrow(data)/1000 >= 1001) {
      print("Multiple predictions not yet optimized computationally... it may time some time to run")
    }
    outriggerfitted <- numeric(ntest)
    for (i in seq_len(ntest)) {
      outriggerfitted[i] <- outrigger_onepoint(formula, data, xtest[i], degree, bandwidth, kernel, lambda, folds, scoregranality, score_df)$prediction
    }
    return(data.frame(x=xtest, y=outriggerfitted))
  }
}

outrigger_onepoint <- function(formula, data, x, degree=0, bandwidth=NULL, kernel="epan", lambda=NULL, folds=5, scoregranality=1, score_df=10){
  #' Performs outrigger regression at a point x
  #'
  #' @param formula formula (as in e.g. \code{lm}).
  #' @param data data frame.
  #' @param x datapoint for local estimation.
  #' @param degree degree of polynomial.
  #' @param bandwidth bandwidth parameter. If NULL a rule-of-thumb bandwidth is selected by least-squares cross-validation for the standard local polynomial estimator.
  #' @param kernel kernel for local polyomial. Default is epanechnikov kernel.
  #' @param lambda orthogonalisation parameter. By default takes \eqn{\lambda = 4(1 + 0.5 \log(n))}.
  #' @param folds number of folds for cross-fitting. Default is 5.
  #' @param score_df number of degrees of freedom for score matching splines used for conditional score estimation. Input can be a numeric value (e.g. 10 on larger datasets, 6 on smaller datasets) or "cross-validation" or "cross-validation-quick", in which case score-matching-CV will be performed (may be computationally costly).
  #' @param scoregranality the number of bins to split covariate-space \eqn{|X-x|\leq\lambda h} for conditional score estimation.
  #'
  #' @return A list containing:
  #'   \describe{
  #'     \item{\code{prediction}}{The outrigger local polynomial estimator \eqn{\hat{f}(x)} at \eqn{x}.}
  #'     \item{\code{fitted_vector}}{For the local linear outrigger estimator, the full fitted vector at \eqn{x}}
  #'     \item{\code{standardlocpol_fitted_vector}}{The fitted vector for the standard local polynomial}
  #'   }
  #' @importFrom np npregbw
  #' @importFrom stats quantile
  #' @export
  #'
  #' @noRd

  n <- nrow(data)

  if (is.null(bandwidth)) {
    # rule of thumb - not necessarily optimal
    print("Calculating an automatic 'rule of thumb' bandwidth calculation...")
    print("NB: Bandwidth selection based on least-squares cross validation for standard local polynomial (see np::npregbw)")
    if (degree==0) {
      bandwidth <- npregbw(xdat=data$X, ydat=data$Y, regtype="lc", ckertype="epanechnikov")$bandwidth$x*sqrt(5)
    } else if (degree==1) {
      bandwidth <- npregbw(xdat=data$X, ydat=data$Y, regtype="ll", ckertype="epanechnikov")$bandwidth$x*sqrt(5)
    } else {
      stop("Degree must be either 0 (local constant) or 1 (local linear)")
    }
    print(paste0("Bandwidth calculation complete: ", bandwidth))
  }
  h <- bandwidth

  if (is.null(lambda)) {
    # rule of thumb - not necessarily optimal
    lambda <- 4 + log(n, base=10)
  }

  if (is.null(scoregranality)) {
    scoregranality <- max(floor(lambda/4),1)
  }
  M <- scoregranality
  if (score_df=="cross-validation-quick") score_df_vec=numeric(M)

  #if (!is.numeric()) {
    #print("NB: Score-matching-cross-validation chosen for fitting at every fold/bin level. This may be rather computationally costly... a user-chosen score_df (and using score_plots) will speed up computation.")
  #}

  if (length(all.vars(formula))==2) {
    response_var <- as.character(formula[[2]])
    covariate_var <- as.character(formula[[3]])
    colnames(data) <- ifelse(colnames(data) == response_var, "Y",
                             ifelse(colnames(data) == covariate_var, "X", colnames(data)))
  } else {
    stop("Currently only one-dimensional regression supported.")
  }

  kernel_fnc <- switch(kernel,
         "epan" = function(u) 1/h*3/4*pmax(0,1-((u-x)/h)^2),
         "gaussian" = function(u) 1/(sqrt(2*pi)*h)*exp(-1/2*((u-x)/h)^2),
         "rectangular" = function(u) 1/(2*h)*as.numeric(abs(u-x)<h),
         "triangular" = function(u) 1/h*pmax(0,1-abs((u-x)/h)),
         "biweight" = function(u) 1/h*15/16*(1-((u-x)/h)^2)^2*as.numeric(abs(u-x)<h),
         stop("Unknown kernel: ", kernel)
         )

  if (degree==0) {
    locpol_predict <- nadaraya_watson_fitted_predict
    locpol_residuals <- nadaraya_watson_residuals
  } else if (degree==1) {
    locpol_predict <- loclin_fitted_predict
    locpol_residuals <- loclin_residuals
  } else {
    stop("degree must be 0 (local constant) or 1 (local linear)")
  }

  data_loc <- data[abs(data$X-x)<=lambda*h,]

  shuffleid <- 1:nrow(data_loc)
  foldid <- cut(shuffleid, breaks = folds, labels = FALSE)

  score_est_k <- list()
  varphi_all <- ind_h_all <- resids_base_all <- numeric(nrow(data_loc))
  for (k in seq_len(folds)) {
    ind_foldk <- foldid==k
    data1 <- data_loc[ind_foldk,]
    data2 <- data_loc[!ind_foldk,]

    K2 <- sapply(data2$X, kernel_fnc)
    dat2_Lambda <- as.numeric( abs(data2$X-x)>h )
    mu0 <- mean(K2)/mean(dat2_Lambda)

    dat1fitted <- locpol_predict(data2$X, data2$Y, h, data1$X, data1$Y)$fitted
    dat1resids <- data1$Y - dat1fitted
    dat1_Lambda <- as.numeric( abs(data1$X-x)>h )
    dat1_h <- (abs(data1$X-x)<=h)
    c0 <- mean(dat1_Lambda*dat1resids)/mean(dat1_Lambda)

    K0 <- sapply(data1$X, kernel_fnc)

    varphi <- K0 - mu0*dat1_Lambda

    # Score function estimation
    data2resids <- locpol_residuals(data2$X, data2$Y, h)
    data2_broken_out <- data.frame(resids=data2resids, X=data2$X)
    #res <- build_subsets_by_quantiles(data2_broken_out, x, h, lambda, M)
    res <- build_subsets_by_quantiles(data2_broken_out, M)
    if (is.numeric(score_df)) {
      for (m in seq_len(M)) res$score[[m]] <- spline_score_RCPP(res$subsets[[m]]$resids, df = score_df)
    } else if (score_df=="cross-validation-quick") {
      for (m in seq_len(M)) {
        if (k==1) score_df_vec[m] = cv_spline_score(res$subsets[[m]]$resids, df=2:10, nfolds=5L)$df_1se
        res$score[[m]] <- spline_score_RCPP(res$subsets[[m]]$resids, df = score_df_vec[m])
      }
    } else {
      for (m in seq_len(M)) {
        score_df_val = cv_spline_score(res$subsets[[m]]$resids, df=2:10, nfolds=5L)$df_1se
        res$score[[m]] <- spline_score_RCPP(res$subsets[[m]]$resids, df = score_df_val)
      }
    }
    score_est_k[[k]] <- res

    varphi_all[ind_foldk] <- varphi
    ind_h_all[ind_foldk] <- dat1_h # T/F
    resids_base_all[ind_foldk] <- dat1resids - c0*dat1_Lambda

    score_ind <- sapply(seq_len(M), function(m) {
      (data1$X >= res$breaks[m]) & (data1$X < res$breaks[m + 1])
    })

  }

  indicators <- sapply(seq_len(M), function(m) {
    (data_loc$X >= res$breaks[m]) & (data_loc$X < res$breaks[m + 1])
  })

  f0 <- locpol_predict(data$X, data$Y, h, x, NA)$fullfitted

    gam0 <- f0
    gam0old <- Inf; acheck <- 0
    while ((abs(gam0[1]-gam0old)>1e-5) & (acheck<=100)) {
      acheck <- acheck + 1
      gam0old <- gam0

      rho_1n <- numeric(nrow(data_loc))
      for (k in seq_len(folds)) {
        ind_foldk <- foldid==k
        ind_foldk_hin <- ind_foldk & ind_h_all==1
        ind_foldk_hout <- ind_foldk & ind_h_all==0
        for (m in 1:M) {
            rho_1n[ind_foldk_hout & indicators[,m]] <- score_est_k[[k]]$score[[m]]$rho(resids_base_all[ind_foldk_hout & indicators[,m]])
            rho_1n[ind_foldk_hin & indicators[,m]] <- score_est_k[[k]]$score[[m]]$rho(data_loc$Y[ind_foldk_hin & indicators[,m]] - gam0)
        }
      }
      num <- varphi_all * rho_1n
      den <- varphi_all * ind_h_all*rho_1n^2
      gam0 <- gam0 - sum(num)/sum(den)

      }
      if (acheck==101) gam0 <- (gam0+gam0old)/2

  output <- list(prediction=gam0[1], fitted_vector=gam0, standardlocpol_fitted_vector=f0,
                 score_plot_metadata = list(all_scores=score_est_k ,
                                             all_pilot_resids=resids_base_all),
                 bwd = h, lambda=lambda
                 )

  print(paste0("Predicted value at x=", x,": ", output$prediction))

  return(output)
}

outrigger_splinescoreinitalizer <- function(formula, data, x, degree=0, bandwidth=NULL, kernel="epan", lambda=NULL, folds=5, scoregranality=1, score_df=10){
  #' Performs outrigger regression at a point x. Splines used for global regression function intializer.
  #'
  #' @param formula formula (as in e.g. \code{lm}).
  #' @param data data frame.
  #' @param x datapoint for local estimation.
  #' @param degree degree of polynomial.
  #' @param bandwidth bandwidth parameter. If NULL a rule-of-thumb bandwidth is selected by least-squares cross-validation for the standard local polynomial estimator.
  #' @param kernel kernel for local polyomial. Default is epanechnikov kernel.
  #' @param lambda orthogonalisation parameter. By default takes \eqn{\lambda = 4(1 + 0.5 \log(n))}.
  #' @param folds number of folds for cross-fitting. Default is 5.
  #' @param score_df number of degrees of freedom for score matching splines used for conditional score estimation. Input can be a numeric value (e.g. 10 on larger datasets, 6 on smaller datasets) or "cross-validation" or "cross-validation-quick", in which case score-matching-CV will be performed (may be computationally costly).
  #' @param scoregranality the number of bins to split covariate-space \eqn{|X-x|\leq\lambda h} for conditional score estimation.
  #'
  #' @return A list containing:
  #'   \describe{
  #'     \item{\code{prediction}}{The outrigger local polynomial estimator \eqn{\hat{f}(x)} at \eqn{x}.}
  #'     \item{\code{fitted_vector}}{For the local linear outrigger estimator, the full fitted vector at \eqn{x}}
  #'     \item{\code{standardlocpol_fitted_vector}}{The fitted vector for the standard local polynomial}
  #'   }
  #' @importFrom np npregbw
  #' @importFrom mgcv gam
  #' @importFrom stats quantile
  #' @export
  #'
  #' @noRd

  n <- nrow(data)


  if (is.null(bandwidth)) {
    # rule of thumb - not necessarily optimal
    print("Calculating an automatic 'rule of thumb' bandwidth calculation...")
    print("NB: Bandwidth selection based on least-squares cross validation for standard local polynomial (see np::npregbw)")
    if (degree==0) {
      bandwidth <- npregbw(xdat=data$X, ydat=data$Y, regtype="lc", ckertype="epanechnikov")$bandwidth$x*sqrt(5)
    } else if (degree==1) {
      bandwidth <- npregbw(xdat=data$X, ydat=data$Y, regtype="ll", ckertype="epanechnikov")$bandwidth$x*sqrt(5)
    } else {
      stop("Degree must be either 0 (local constant) or 1 (local linear)")
    }
    print(paste0("Bandwidth calculation complete: ", bandwidth))
  }
  h <- bandwidth

  if (is.null(lambda)) {
    # rule of thumb - not necessarily optimal
    lambda <- 4 + log(n, base=10)
  }

  if (is.null(scoregranality)) {
    scoregranality <- max(floor(lambda/4),1)
  }
  M <- scoregranality
  if (score_df=="cross-validation-quick") score_df_vec=numeric(M)

  #if (!is.numeric()) {
  #print("NB: Score-matching-cross-validation chosen for fitting at every fold/bin level. This may be rather computationally costly... a user-chosen score_df (and using score_plots) will speed up computation.")
  #}

  if (length(all.vars(formula))==2) {
    response_var <- as.character(formula[[2]])
    covariate_var <- as.character(formula[[3]])
    colnames(data) <- ifelse(colnames(data) == response_var, "Y",
                             ifelse(colnames(data) == covariate_var, "X", colnames(data)))
  } else {
    stop("Currently only one-dimensional regression supported.")
  }

  kernel_fnc <- switch(kernel,
                       "epan" = function(u) 1/h*3/4*pmax(0,1-((u-x)/h)^2),
                       "gaussian" = function(u) 1/(sqrt(2*pi)*h)*exp(-1/2*((u-x)/h)^2),
                       "rectangular" = function(u) 1/(2*h)*as.numeric(abs(u-x)<h),
                       "triangular" = function(u) 1/h*pmax(0,1-abs((u-x)/h)),
                       "biweight" = function(u) 1/h*15/16*(1-((u-x)/h)^2)^2*as.numeric(abs(u-x)<h),
                       stop("Unknown kernel: ", kernel)
  )

  if (degree==0) {
    locpol_predict <- nadaraya_watson_fitted_predict
    locpol_residuals <- nadaraya_watson_residuals
  } else if (degree==1) {
    locpol_predict <- loclin_fitted_predict
    locpol_residuals <- loclin_residuals
  } else {
    stop("degree must be 0 (local constant) or 1 (local linear)")
  }

  data_loc <- data[abs(data$X-x)<=lambda*h,]

  shuffleid <- 1:nrow(data_loc)
  foldid <- cut(shuffleid, breaks = folds, labels = FALSE)

  score_est_k <- list()
  varphi_all <- ind_h_all <- resids_base_all <- numeric(nrow(data_loc))
  for (k in seq_len(folds)) {
    ind_foldk <- foldid==k
    data1 <- data_loc[ind_foldk,]
    data2 <- data_loc[!ind_foldk,]

    K2 <- sapply(data2$X, kernel_fnc)
    dat2_Lambda <- as.numeric( abs(data2$X-x)>h )
    mu0 <- mean(K2)/mean(dat2_Lambda)

    dat1fitted <- predict(gam(Y~s(X,bs="cr"), data=data2), data1)
    dat1resids <- data1$Y - dat1fitted
    dat1_Lambda <- as.numeric( abs(data1$X-x)>h )
    dat1_h <- (abs(data1$X-x)<=h)
    c0 <- mean(dat1_Lambda*dat1resids)/mean(dat1_Lambda)

    K0 <- sapply(data1$X, kernel_fnc)

    varphi <- K0 - mu0*dat1_Lambda

    # Score function estimation
    data2resids <- gam(Y~s(X,bs="cr"), data=data2)$residuals  #locpol_residuals(data2$X, data2$Y, h) #0.15
    data2_broken_out <- data.frame(resids=data2resids, X=data2$X)
    res <- build_subsets_by_quantiles(data2_broken_out, M)
    if (is.numeric(score_df)) {
      for (m in seq_len(M)) res$score[[m]] <- spline_score_RCPP(res$subsets[[m]]$resids, df = score_df)
    } else if (score_df=="cross-validation-quick") {
      for (m in seq_len(M)) {
        if (k==1) score_df_vec[m] = cv_spline_score(res$subsets[[m]]$resids, df=2:10, nfolds=5L)$df_1se
        res$score[[m]] <- spline_score_RCPP(res$subsets[[m]]$resids, df = score_df_vec[m])
      }
    } else {
      for (m in seq_len(M)) {
        score_df_val = cv_spline_score(res$subsets[[m]]$resids, df=2:10, nfolds=5L)$df_1se
        res$score[[m]] <- spline_score_RCPP(res$subsets[[m]]$resids, df = score_df_val)
      }
    }
    score_est_k[[k]] <- res

    varphi_all[ind_foldk] <- varphi
    ind_h_all[ind_foldk] <- dat1_h # T/F
    resids_base_all[ind_foldk] <- dat1resids - c0*dat1_Lambda

    score_ind <- sapply(seq_len(M), function(m) {
      (data1$X >= res$breaks[m]) & (data1$X < res$breaks[m + 1])
    })

  }

  indicators <- sapply(seq_len(M), function(m) {
    (data_loc$X >= res$breaks[m]) & (data_loc$X < res$breaks[m + 1])
  })

  f0 <- locpol_predict(data$X, data$Y, h, x, NA)$fullfitted

  gam0 <- f0
  gam0old <- Inf; acheck <- 0
  while ((abs(gam0[1]-gam0old)>1e-5) & (acheck<=100)) {
    acheck <- acheck + 1
    gam0old <- gam0

    rho_1n <- numeric(nrow(data_loc))
    for (k in seq_len(folds)) {
      ind_foldk <- foldid==k
      ind_foldk_hin <- ind_foldk & ind_h_all==1
      ind_foldk_hout <- ind_foldk & ind_h_all==0
      for (m in 1:M) {
        rho_1n[ind_foldk_hout & indicators[,m]] <- score_est_k[[k]]$score[[m]]$rho(resids_base_all[ind_foldk_hout & indicators[,m]])
        rho_1n[ind_foldk_hin & indicators[,m]] <- score_est_k[[k]]$score[[m]]$rho(data_loc$Y[ind_foldk_hin & indicators[,m]] - gam0)
      }
    }
    num <- varphi_all * rho_1n
    den <- varphi_all * ind_h_all*rho_1n^2
    gam0 <- gam0 - sum(num)/sum(den)

  }
  if (acheck==101) gam0 <- (gam0+gam0old)/2

  output <- list(prediction=gam0[1], fitted_vector=gam0, standardlocpol_fitted_vector=f0,
                 score_plot_metadata = list(all_scores=score_est_k ,
                                            all_pilot_resids=resids_base_all),
                 bwd = h, lambda=lambda
  )

  return(output)
}

build_subsets_by_quantiles <- function(data, M) {
  # Ensure M is at least 1
  stopifnot(M >= 1)

  # Calculate the quantiles of data$X
  quantiles <- quantile(data$X, probs = seq(0, 1, length.out = M + 1))

  # Generate subsets based on these quantiles
  subsets <- lapply(seq_len(M), function(i) {
    data[data$X >= quantiles[i] & data$X < quantiles[i + 1], ]
  })

  # Return the subsets and the break points (quantiles)
  return(list(subsets = subsets, breaks = quantiles))
}
