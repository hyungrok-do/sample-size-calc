// A Super Learner for the nuisance models of AIPW and TMLE: a cross-validated
// ensemble that needs the adjustment variables but not their functional form.
#pragma once
#include "model.hpp"
#include <string>
#include <vector>

namespace ssc::dag {
// The library, in the order of every weight vector: the training mean, a
// main-terms GLM, a GLM with all pairwise products, additive restricted cubic
// splines (knots at covariate quantiles), gradient-boosted depth-2 trees with
// early stopping, and kernel ridge regression on random Fourier features.
extern const std::vector<std::string> LEARNERS;
struct Learned {
  Vec prediction;  // for every row of the data
  Vec weights;     // per learner, summing to one; zero for a learner that failed
};
// Fits to `rows` and predicts every row. Each learner's 5-fold cross-validated
// predictions are combined by non-negative least squares, normalized to a
// convex combination, and the learners with weight are refitted on all `rows`.
// family "bernoulli" uses logistic GLMs and logistic boosting and returns
// probabilities; any other family uses squared-error learners.
Learned super_learner(const Data &data, const std::vector<int> &variables,
                      const Vec &response, const std::vector<int> &rows,
                      const std::string &family);
} // namespace ssc::dag
