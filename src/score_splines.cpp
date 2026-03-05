#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector ng_pseudo_response_rcpp(NumericVector x, NumericVector w = NumericVector::create()) {
  int n = x.size();

  // Default weights = 1
  if (w.size() == 0) {
    w = NumericVector(n, 1.0);
  }

  // Order x
  IntegerVector order0 = seq(0, n - 1);
  std::sort(order0.begin(), order0.end(), [&](int i, int j) { return x[i] < x[j]; });

  // Sort x and w
  NumericVector xs(n), ws(n);
  for (int i = 0; i < n; i++) {
    xs[i] = x[order0[i]];
    ws[i] = w[order0[i]];
  }

  // Interval widths
  NumericVector h(n - 1);
  for (int i = 0; i < n - 1; i++) {
    h[i] = xs[i + 1] - xs[i];
  }

  // wih and wh
  NumericVector wih(n - 1), wh(n - 1);
  for (int i = 0; i < n - 2; i++) {
    wih[i] = ws[i] / h[i];
    wh[i] = ws[i] * h[i];
  }
  wih[n - 2] = (ws[n - 2] + ws[n - 1]) / h[n - 2];
  wh[n - 2] = (ws[n - 2] - ws[n - 1] / 2.0) * h[n - 2];

  // a_vec
  vec a_vec(n, fill::zeros);
  for (int i = 0; i < n - 1; i++) {
    a_vec[i]     += wih[i];
    a_vec[i + 1] -= wih[i];
  }

  // c_vec
  vec c_vec(n - 2);
  for (int i = 0; i < n - 2; i++) {
    c_vec[i] = (wh[i] + 2.0 * wh[i + 1]) / 3.0;
  }

  // ih
  NumericVector ih(n - 1);
  for (int i = 0; i < n - 1; i++) {
    ih[i] = 1.0 / h[i];
  }

  // R matrix
  mat R(n - 2, n - 2, fill::zeros);
  for (int i = 0; i < n - 2; i++) {
    R(i, i) = 2.0 * (h[i] + h[i + 1]) / 3.0;
    if (i > 0) R(i, i - 1) = h[i] / 3.0;
    if (i < n - 3) R(i, i + 1) = h[i + 1] / 3.0;
  }

  // Q matrix
  mat Q(n, n - 2, fill::zeros);
  for (int j = 0; j < n - 2; j++) {
    Q(j, j)       = ih[j];
    Q(j + 1, j)   = -(ih[j] + ih[j + 1]);
    Q(j + 2, j)   = ih[j + 1];
  }

  // Solve R z = c_vec
  vec z = solve(R, c_vec);

  // Compute y
  vec y_sorted = a_vec + Q * z;
  for (int i = 0; i < n; i++) {
    y_sorted[i] /= ws[i];
  }

  // Restore original order
  NumericVector y(n);
  for (int i = 0; i < n; i++) {
    y[order0[i]] = y_sorted[i];
  }

  return y;
}

