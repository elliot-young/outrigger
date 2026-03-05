nadaraya_watson_residuals <- function(X, Y, h) {
  # Check input lengths
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length.")
  }

  n <- length(X)
  m_hat <- numeric(n)  # Predicted values

  # Epanechnikov kernel function
  epanechnikov_kernel <- function(u) {
    k <- 0.75 * (1 - u^2)
    k[abs(u) > 1] <- 0
    return(k)
  }

  # Compute Nadaraya-Watson estimator at each point
  for (i in 1:n) {
    u <- (X[i] - X) / h
    weights <- epanechnikov_kernel(u)
    if (sum(weights) == 0) {
      m_hat[i] <- NA  # Avoid division by zero
    } else {
      m_hat[i] <- sum(weights * Y) / sum(weights)
    }
  }

  # Residuals
  residuals <- Y - m_hat
  return(residuals)
}

nadaraya_watson_fitted_predict <- function(X, Y, h, Xnew, Ynew) {
  # Check input lengths
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length.")
  }

  n <- length(X)

  nnew <- length(Xnew)
  m_hat <- numeric(nnew)

  # Epanechnikov kernel function
  epanechnikov_kernel <- function(u) {
    k <- 0.75 * (1 - u^2)
    k[abs(u) > 1] <- 0
    return(k)
  }

  # Compute Nadaraya-Watson estimator at each point
  for (i in 1:nnew) {
    u <- (Xnew[i] - X) / h
    weights <- epanechnikov_kernel(u)
    if (sum(weights) == 0) {
      m_hat[i] <- NA  # Avoid division by zero
    } else {
      m_hat[i] <- sum(weights * Y) / sum(weights)
    }
  }

  # Residuals
  fitted <- m_hat
  return(return(list(fitted=fitted, fullfitted=fitted)))
}

loclin_residuals <- function(X, Y, h) {
  # Check input lengths
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length.")
  }

  n <- length(X)
  m_hat <- numeric(n)  # Predicted values

  # Epanechnikov kernel function
  epanechnikov_kernel <- function(u) {
    k <- 0.75 * (1 - u^2)
    k[abs(u) > 1] <- 0
    return(k)
  }

  # Compute Nadaraya-Watson estimator at each point
  for (i in 1:n) {
    u <- (X[i] - X) / h
    weights <- epanechnikov_kernel(u)
    if (sum(weights) == 0) {
      m_hat[i] <- NA  # Avoid division by zero
    } else {
      S0 <- sum(weights)
      S1 <- sum(weights * u)
      S2 <- sum(weights * u^2)
      Y0 <- sum(weights * Y)
      Y1 <- sum(weights * u * Y)
      m_hat[i] <- (S2*Y0-S1*Y1)/(S0*S2-S1^2)
    }
  }

  # Residuals
  residuals <- Y - m_hat
  return(residuals)
}

loclin_fitted_predict <- function(X, Y, h, Xnew, Ynew) {
  # Check input lengths
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length.")
  }

  n <- length(X)

  nnew <- length(Xnew)
  m_hat <- array(NA, dim=c(nnew,2)) #numeric(nnew)

  # Epanechnikov kernel function
  epanechnikov_kernel <- function(u) {
    k <- 0.75 * (1 - u^2)
    k[abs(u) > 1] <- 0
    return(k)
  }

  # Compute Nadaraya-Watson estimator at each point
  for (i in 1:nnew) {
    u <- (Xnew[i] - X) / h
    weights <- epanechnikov_kernel(u)
    if (sum(weights) == 0) {
      m_hat[i] <- NA  # Avoid division by zero
    } else {
      S0 <- sum(weights)
      S1 <- sum(weights * u)
      S2 <- sum(weights * u^2)
      Y0 <- sum(weights * Y)
      Y1 <- sum(weights * u * Y)
      m_hat[i,1] <- (S2*Y0-S1*Y1)/(S0*S2-S1^2)
      m_hat[i,2] <- (S0*Y1-S1*Y0)/(S0*S2-S1^2)
    }
  }

  # Residuals
  fitted <- m_hat
  return(list(fitted=fitted[,1], fullfitted=fitted))
}
