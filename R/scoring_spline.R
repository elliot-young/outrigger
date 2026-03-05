spline_score_RCPP <- function(x, df=5, tol=1e-3, nmax=NULL){
  #' Univariate score estimation via the smoothing spline method of Cox 1985 and Ng 1994.
  #'
  #' @param x vector of datapoints
  #' @param df vector of smoothing parameters for the
  #'     non-parametric score estimator, corresponding to the
  #'     effective degrees of freedom for a smoothing spline.
  #' @param tol numeric tolerance, minimum distance between neighbouring points,
  #'     to avoid singularities.
  #' @param nmax if specified, overrides tol as maximal number of unique points.
  #'
  #' @return score function "rho" and derivative "drho", which take vector
  #' input and yield a vector of score estimates corresponding to each df (in a
  #' list if there are multiple df values). Also output the vector "df".
  #' @noRd

  if (!requireNamespace("graphics", quietly = TRUE)) {
    warning("Package \"graphics\" required for sort_bin. Matrices may be
            singular.")
    w <- rep(1,length(x))
    x_sort <- x
  } else {
    bin <- sort_bin(x=x, tol=tol, nmax=nmax)
    x_sort <- bin$x_sort
    w <- bin$w

    if (length(x_sort) < 4){
      warning("smooth.spline requires at least 4 unique x values.
                Binning would violate this, so not binning.")
      w <- rep(1,length(x))
      x_sort <- x
    }
  }

  pseudo_y <- ng_pseudo_response_rcpp(x=x_sort, w=w)

  rho <- function(x){
    out <- lapply(df, function(dfval){
      sm <- stats::smooth.spline(x=x_sort, y=pseudo_y, w=w, df=dfval)
      return(as.vector(stats::predict(sm, x)$y))
    })
    # if only one df specified, unlist
    if (length(out)==1){
      out <- out[[1]]
    }
    return(out)
  }

  drho <- function(x){
    out <- lapply(df, function(dfval){
      sm <- stats::smooth.spline(x=x_sort, y=pseudo_y, w=w, df=dfval)
      return(as.vector(stats::predict(sm, x, deriv=1)$y))
    })
    # if only one df specified, unlist
    if (length(out)==1){
      out <- out[[1]]
    }
    return(out)
  }

  return(list("rho" = rho,
              "drho" = drho,
              "df" = df))
}

cv_spline_score <- function(x, df=2:10, nfolds=5L, tol=1e-3, nmax=NULL){
  #' Cross-validation for spline_score.
  #'
  #' @param x vector of errors (residuals)
  #' @param df vector of smoothing parameters for the
  #'     score estimator, corresponding to the
  #'     effective degrees of freedom for a smoothing spline.
  #' @param nfolds integer number of cross-validation folds.
  #' @param tol numeric tolerance, minimum distance between neighbouring points,
  #'     to avoid singularities.
  #' @param nmax if specified, overrides tol as maximal number of unique points.
  #'
  #' @return list of 5 elements: df vector, cv vector of corresponding
  #'     cross-validation scores, se vector of standard error estimates,
  #'     df_min cross-validation minimiser, df_1se largest smoothing
  #'     parameter within CV score within one standard error of df_min.
  #' @noRd

  n <- length(x)
  ndf <- length(df)

  foldid <- sample(rep(seq(nfolds), length.out=n))
  cv_folds <- matrix(NA,nrow=nfolds,ncol=ndf)

  for (fold in seq(nfolds)){

    which <- (foldid == fold)
    spline <- spline_score_RCPP(x=x[!which], df=df, tol=tol, nmax=nmax)

    # Evaluate score estimate on holdout fold
    score <- spline$rho(x[which])
    dscore <- spline$drho(x[which])
    # Compute CV scores for each df
    cv_folds[fold, ] <- mapply(function(x,y){mean(x^2+2*y)},
                               x = score,
                               y = dscore)
  }

  cv <- colMeans(cv_folds)
  se <- apply(cv_folds,2,stats::sd)/sqrt(nfolds)

  dfmin_index <- which.min(cv)
  cv_target <- cv[dfmin_index]+se[dfmin_index]
  df1se_index <- which(df==min(df[cv < cv_target]))

  df_min <- df[dfmin_index]
  df_1se <- df[df1se_index]

  return(list('df'=df,
              'cv'=cv,
              'se'=se,
              'df_min'=df_min,
              'df_1se'=df_1se))

}

sort_bin <- function(x, tol=1e-5, nmax=NULL){
  #' Sort and bin x within a specified tolerance, using hist().
  #'
  #' @param x vector of covariates.
  #' @param tol numeric tolerance, minimum distance between neighbouring points,
  #'     to avoid singularities.
  #' @param nmax if specified, overrides tol as maximal number of unique points.
  #'
  #' @return list with three elements. x_sort is sorted and binned x,
  #'      w is a vector of weights corresponding to the frequency of each bin,
  #'      order0 is a vector specifying the order0ing of x into the binned values
  #'      sort_x.
  #' @noRd


  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Package \"graphics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }


  if (is.null(nmax)){
    br <- ceiling((max(x)-min(x))/tol) # number of bins
  } else{
    br <- nmax
  }

  br <- min(br, 1e+6) # hist caps br to 1e+6 anyway, so this just avoid a warning message

  hist <- graphics::hist(x, br, right=FALSE, plot=FALSE) # assign elements of x to bins
  counts <- hist$counts
  mids <- hist$mids
  breaks <- hist$breaks

  # remove empty bins
  w <- counts[counts>0] # frequencies in non-empty bins
  x_sort <- mids[counts>0] # midpoints of non-empty bins
  pos_breaks <- breaks[c(1,1+which(counts>0))] # breakpoints for non-empty bins

  order0 <- findInterval(x, pos_breaks)

  return(list("x_sort"=x_sort, "w"=w, "order0"=order0))

}

