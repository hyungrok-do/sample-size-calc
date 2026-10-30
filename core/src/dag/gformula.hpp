// The parametric g-formula through the DAG: the analysis used when no simpler
// test identifies a marginal target.
#pragma once
#include "model.hpp"

namespace ssc::dag {
// A variable the analysis refits, with the terms of its conditional model.
struct ChainModel {
  int node = 0;
  std::vector<Term> terms;
};
struct GFormula {
  double estimate = 0, testEstimate = 0, se = 0;
};
// `chain` lists, in generation order, the variables the comparison variable
// affects on the way to the outcome, ending with the outcome. Each is refitted
// to the study's data with its own family (censored event times by Weibull
// maximum likelihood); every other variable keeps its observed values. The
// target is the average, over the study's rows, of the outcome functional
// (mean, P(Y = 1), S(horizon) or RMST to the horizon) with the comparison
// variable set to each level and the chain drawn from the fitted models with
// common random numbers. Its standard error adds the rows' variation to the
// fitted parameters' (a delta method evaluated at +/- one Cholesky column of
// each model's covariance).
GFormula gformula(const Graph &graph, const Data &data,
                  const std::vector<ChainModel> &chain, int exposure,
                  double reference, double active, const std::string &measure,
                  double horizon);
} // namespace ssc::dag
