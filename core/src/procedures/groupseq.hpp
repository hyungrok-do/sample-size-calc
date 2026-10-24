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

/// Error spent between information fractions s and t, 0 ≤ s < t ≤ 1 (the whole level by t = 1).
using SpendFn = std::function<double(double s, double t)>;

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
  SpendFn spend;     // alpha spent between two information fractions (Spending)
  double delta = 0;  // b_k = C·t_k^(delta − ½) (WangTsiatis): 0 = O'Brien–Fleming, ½ = Pocock
};

struct Futility {
  enum class Kind { None, BetaSpending, Fixed };
  Kind kind = Kind::None;
  SpendFn spend;          // beta spent between two information fractions (BetaSpending)
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
  int grid = 32;            // r of the grid its probabilities are computed on (see grid_resolution)
  bool has_futility() const { return !a.empty(); }
};

/// The grid's r for a design at level alpha with information fractions t: the smallest of 32 (as Jennison &
/// Turnbull advise), 48, 64, 96 and 128 that keeps the last boundary within about 2e-5 of a grid with r = 256 for up
/// to 20 equally spaced analyses: 32 down to an alpha of 1e-4, 48 to 1e-8, 64 to 1e-12, 96 to 1e-20 and 128 below.
/// The smaller the alpha, the farther out the boundaries lie, where the sparse tails of a coarse grid lose
/// accuracy: with 32 and up to 20 analyses the error was 1e-4 at 1e-7 and 3e-4 at 1e-12, N was 2 more at 1e-5
/// (Wang–Tsiatis, 20 analyses), and at 1e-20 the last nominal level exceeded alpha.  Two analyses closer than 20
/// equally spaced ones make the step from one to the next narrower than that: below 1e-8 such designs take 128.
/// An alpha-spending function that spends less than a thousandth of alpha at the last analysis takes the r of a
/// thousand times that (with ρ = 1e-30 and 20 analyses, N was 2 short with 32), and so does a β-spending function
/// of futility boundaries that leaves less than a thousandth of beta for it (with a β-spending ρ of 1e-30 and 20
/// analyses, N was 2 or 4 off).  The work grows with r².
int grid_resolution(double alpha, const std::vector<double>& t, const Efficacy& eff, const Futility& fut = {});

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
/// steps and may throw to cancel.  r = 0 takes grid_resolution(alpha, t, eff).
Design make_design(const std::vector<double>& t, double alpha, bool two_sided, const Efficacy& eff,
                   const Futility& fut, const std::function<void()>& check = {}, int r = 0);

/// Exit probabilities under drift eta.  `scale` multiplies every boundary (a
/// statistic whose standard deviation under H1 differs from 1 can be handled by
/// scaling its boundaries by 1/SD and its drift by 1/SD).  r = 0 takes the design's grid.
ExitProbs exit_probs(const Design& d, double eta, double scale = 1, int r = 0);

/// Drift giving the requested overall power (rejection in the direction of the
/// drift); boundaries are multiplied by `scale` as in exit_probs.
double drift_for_power(const Design& d, double power, double scale = 1, int r = 0);

}  // namespace ssc::gs
