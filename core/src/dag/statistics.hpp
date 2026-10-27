#pragma once
#include "model.hpp"
#include <cmath>

namespace ssc::dag {
struct Moments {
  int n = 0;
  double mean = 0, variance = 0;
};
Moments moments(const Vec &values);
struct Matrix {
  Vec x, scales;
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
Fit fit_cox(const Matrix &matrix, const Vec &times, const Vec &events,
            int target);
struct KM {
  double survival = 1, survivalVariance = 0, rmst = 0, rmstVariance = 0;
};
KM kaplan_meier(const Vec &times, const Vec &events, double horizon);
double p_value(double statistic, double df, const std::string &side);
} // namespace ssc::dag
