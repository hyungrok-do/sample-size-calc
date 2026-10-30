#pragma once
#include "model.hpp"
#include <cmath>

namespace ssc::dag {
struct Moments {
  int n = 0;
  double mean = 0, variance = 0;
};
Moments moments(const Vec &values);
// Columns are centred and scaled (the intercept is left at 1): column j holds
// (value - means[j]) / scales[j].
struct Matrix {
  Vec x, scales, means;
  int n = 0, p = 0;
};
Matrix model_matrix(const Data &data, const std::vector<Term> &terms,
                    bool intercept);
struct Fit {
  double coefficient = 0, se = 0, df = INFINITY;
};
Fit fit_linear(const Matrix &matrix, const Vec &response, int target);
Fit fit_glm(const Matrix &matrix, const Vec &response, int target,
            const std::string &family);
// Means for every row of `matrix` from a model fitted to `rows` only: least
// squares ("gaussian") or a canonical/log-link GLM ("bernoulli", "poisson",
// "gamma").
Vec fitted_means(const Matrix &matrix, const Vec &response,
                 const std::vector<int> &rows, const std::string &family);
Fit fit_cox(const Matrix &matrix, const Vec &times, const Vec &events,
            int target);
// A whole fitted model on the matrix's scaled columns: coefficients and their
// covariance (multiplied by `dispersion`, the residual variance of least
// squares or the Pearson dispersion of a Gamma GLM, and 1 otherwise).
struct ModelFit {
  Vec beta, covariance;
  double dispersion = 1;
};
// "gaussian" (least squares), "bernoulli", "poisson" or "gamma" (log link).
ModelFit fit_model(const Matrix &matrix, const Vec &response,
                   const std::string &family);
// Weibull proportional hazards with right censoring, H(t | x) = t^k exp(x'b):
// `beta` holds b and then log k, unless `exponential` fixes k = 1.
ModelFit fit_weibull(const Matrix &matrix, const Vec &times, const Vec &events,
                     bool exponential);
struct KM {
  double survival = 1, survivalVariance = 0, rmst = 0, rmstVariance = 0;
};
KM kaplan_meier(const Vec &times, const Vec &events, double horizon);
double p_value(double statistic, double df, const std::string &side);
} // namespace ssc::dag
