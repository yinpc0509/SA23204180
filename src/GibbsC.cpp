#include <Rcpp.h>
using namespace Rcpp;

// Function to perform Gibbs sampling in C++
// [[Rcpp::export]]
NumericMatrix gibbsSampler(int num_samples, int n, double a, double b) {
  NumericMatrix samples(num_samples, 2);

  // Initial values
  double x = 0;
  double y = 0.5;

  for (int i = 0; i < num_samples; ++i) {
    // Sample x from Binomial(n, y)
    x = R::rbinom(n, y);

    // Sample y from Beta(x + a, n - x + b)
    y = R::rbeta(x + a, n - x + b);

    samples(i, 0) = x;
    samples(i, 1) = y;
  }

  return samples;
}

