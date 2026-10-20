// Group-sequential designs: efficacy boundaries from alpha-spending functions
// or the Wang–Tsiatis family, optional futility boundaries (beta-spending or
// fixed z-values, binding or non-binding), exit probabilities and the drift
// needed for a given power.
//
// Probabilities follow the canonical joint distribution of the z-statistics
// Z_1..Z_K at information fractions t_1 < ... < t_K = 1 with drift eta (the
// mean of Z_K).  They are computed with the recursive numerical integration of
// Armitage, McPherson & Rowe (1969) on the grid of Jennison & Turnbull (2000,
// ch. 19): Simpson's rule on about 12r points per analysis, O(K r²) work.
#pragma once

#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace ssc::gs {

/// Cumulative error spent at information fraction t ∈ [0, 1].
using SpendFn = std::function<double(double t)>;

/// One-sided spending function at total level `level`: "obf" (Lan–DeMets
/// O'Brien–Fleming type), "pocock" (Lan–DeMets Pocock type), "power"
/// (Kim–DeMets level·t^param), "hsd" (Hwang–Shih–DeCani with γ = param).
SpendFn spending(const std::string& name, double level, double param);

/// Alpha-spending function of a design with total level alpha.  A symmetric
/// two-sided design spends alpha/2 per side with the one-sided function, the
/// convention of rpact, gsDesign and ldbounds (it only matters for "obf").
SpendFn alpha_spending(const std::string& name, double alpha, bool two_sided, double param);

struct Efficacy {
  enum class Kind { Spending, WangTsiatis };
  Kind kind = Kind::Spending;
  SpendFn spend;     // cumulative alpha at t (Spending)
  double delta = 0;  // b_k = C·t_k^(delta − ½) (WangTsiatis): 0 = O'Brien–Fleming, ½ = Pocock
};

struct Futility {
  enum class Kind { None, BetaSpending, Fixed };
  Kind kind = Kind::None;
  SpendFn spend;          // cumulative beta at t (BetaSpending)
  double beta = 0.2;      // total type II error (BetaSpending)
  std::vector<double> z;  // z-scale futility boundaries at analyses 1..K−1 (Fixed)
  bool binding = false;   // binding boundaries are used when computing the type I error
};

struct Design {
  std::vector<double> t;            // information fractions, t.back() == 1
  std::vector<double> b;            // efficacy boundaries: reject when Z_k ≥ b_k (|Z_k| ≥ b_k if two-sided)
  std::vector<double> a;            // futility boundaries, stop when Z_k ≤ a_k; empty when none (a.back() == b.back())
  std::vector<double> alpha_spent;  // cumulative type I error spent by the efficacy boundaries
  std::vector<double> beta_spent;   // cumulative type II error spent under eta_design (beta-spending only)
  bool two_sided = false;
  bool binding = false;
  double eta_design = NAN;  // drift under which beta-spending futility boundaries were derived
  bool has_futility() const { return !a.empty(); }
};

struct ExitProbs {
  std::vector<double> up;   // P(stop at analysis k by crossing +b_k)
  std::vector<double> lo;   // P(stop at analysis k by crossing −b_k) (two-sided designs)
  std::vector<double> fut;  // P(stop at analysis k for futility); the last entry is the final non-rejection
  double reject = 0;        // total rejection probability
  double info = 0;          // expected information fraction at stopping (E(N) / N_max)
};

/// Boundaries of a design with total level alpha.  Futility boundaries need a
/// one-sided design; beta-spending futility needs alpha-spending efficacy
/// boundaries and solves for the drift at which the total beta is spent
/// (stored in eta_design).  `check` (optional) is called between expensive
/// steps and may throw to cancel.
Design make_design(const std::vector<double>& t, double alpha, bool two_sided, const Efficacy& eff,
                   const Futility& fut, const std::function<void()>& check = {}, int r = 32);

/// Exit probabilities under drift eta.  `scale` multiplies every boundary (a
/// statistic whose standard deviation under H1 differs from 1 can be handled by
/// scaling its boundaries by 1/SD and its drift by 1/SD).
ExitProbs exit_probs(const Design& d, double eta, double scale = 1, int r = 32);

/// Drift giving the requested overall power (rejection in the direction of the
/// drift); boundaries are multiplied by `scale` as in exit_probs.
double drift_for_power(const Design& d, double power, double scale = 1, int r = 32);

}  // namespace ssc::gs
